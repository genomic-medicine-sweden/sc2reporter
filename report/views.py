from flask_pymongo import PyMongo
from report import app
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, send_from_directory, make_response, escape
from flask_login import LoginManager, UserMixin, login_required, login_user, logout_user, current_user
import json
from bson.objectid import ObjectId
import pprint
from collections import defaultdict
import re
import urllib.parse
from datetime import datetime, timedelta, timezone
import os
import pytz
from skbio import DistanceMatrix
from skbio.tree import nj
import requests
import random
from report.forms import LoginForm
from report.user import User
from operator import itemgetter

login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"

@app.route('/', methods=['GET'])
@login_required
def index():
    show_validation = request.args.get("validation", "")
    search_string = request.args.get("search_string", "")

    query = build_sample_query(search_string=search_string, validation=show_validation)
    samples_cursor = app.config['SAMPLE_COLL'].find(query).sort("_id", -1)
        
    samples = add_significant_variants(samples_cursor)
    lineages_of_concern = app.config["PANGO_LINEAGES_OF_CONCERN"]
    variants_of_significance = set(app.config["VARIANTS_OF_BIOLOGICAL_SIGNIFICANCE"])
    positions_of_significance = set(app.config["POSITIONS_OF_BIOLOGICAL_SIGNIFICANCE"])
    verbosity = request.args.get("verbosity", "simple")
    return render_template('main.html',
                           samples=samples,
                           lineages_of_concern=lineages_of_concern,
                           variants_of_significnace=variants_of_significance,
                           positions_of_significance=positions_of_significance,
                           verbosity=verbosity, search_string=search_string)


@app.route('/reruns')
@login_required
def reruns():
    samples = app.config['SAMPLE_COLL'].find({'validation':{'$ne':1}}).sort("time_added", -1)
    sample_dict = defaultdict(list)
    for s in samples:
        sample_dict[s["sample_id"]].append(s)

    return render_template('reruns.html', samples=sample_dict)

    

@app.route('/dashboard', methods=['GET'])
@login_required
def dashboard():

    criterion_checkboxes = request.args

    success_samples_cursor = app.config['SAMPLE_COLL'].find({'qc.pct_N_bases':{'$lte':app.config["QC_MAX_PCT_N"]}, 'hidden':{'$ne':1}, 'validation':{'$ne':1}}, {'sample_id':1, 'qc':1, 'pangolin':1, 'time_added':1, 'variants':1, 'vcf_filename':1, 'hidden':1, 'collection_date':1, 'selection_criterion':1}).sort("collection_date", -1)
    failed_samples_count = app.config['SAMPLE_COLL'].find({'qc.pct_N_bases':{'$gt':app.config["QC_MAX_PCT_N"]}, 'hidden':{'$ne':1}, 'validation':{'$ne':1}}).count()

    today = datetime.today()

    variant_cursor = app.config['VARIANT_COLL'].find()

    success_samples = list(success_samples_cursor)
    pango_per_date, pango_total_counts = rolling_mean_pango_types(success_samples, criterion_checkboxes)
    pango_count_per_week, weeks = pango_per_week(success_samples, criterion_checkboxes)

    lineages_of_concern = app.config["PANGO_LINEAGES_OF_CONCERN"]
   
    
    return render_template('dashboard.html',
                           samples=success_samples_cursor,
                           failed_samples=failed_samples_count,
                           pango_per_date=pango_per_date,
                           pango_total_counts=pango_total_counts,
                           variants=variant_cursor,
                           pango_count_per_week=pango_count_per_week,
                           weeks=weeks,
                           lineages_of_concern=lineages_of_concern,
                           criterion_checkboxes=criterion_checkboxes
    )


@app.route('/report/<string:sample_id>/<int:max_diff>/<string:verbosity>')
@login_required
def report(sample_id, max_diff, verbosity):
    sample_to_report = app.config['SAMPLE_COLL'].find_one({'_id':ObjectId(sample_id)})

    replicates = []
    if verbosity=="advanced":
        all_samples = app.config['SAMPLE_COLL'].find()
        replicates = app.config['SAMPLE_COLL'].find({'sample_id':sample_to_report["sample_id"]})
    else:
        all_samples = app.config['SAMPLE_COLL'].find({'qc.pct_N_bases':{'$lte':app.config["QC_MAX_PCT_N"]}, 'hidden':{'$ne':1}})

    variants_of_significance = set(app.config["VARIANTS_OF_BIOLOGICAL_SIGNIFICANCE"])
    positions_of_significance = set(app.config["POSITIONS_OF_BIOLOGICAL_SIGNIFICANCE"])    
    tot_samples = all_samples.count()
        

    # Find similar samples (expensive!)
    samples_to_show, variants, report_variants, sample_counts = get_similar_samples(sample_to_report, all_samples, max_diff)
    
    sample_ids = []
    variant_data = {}
    positions = []
    for sample in sorted(samples_to_show, key=lambda k: k["n_diff"]):
        #sid = sample["sample_id"]
        sid = str(sample["_id"]) 
        sample_ids.append(sid)
        variant_data[sid] = {}
        variant_data[sid]["variants"] = {}
        variant_data[sid]["sample_data"] = sample
        for var in variants:
            pos,bases = var.split('_')
            ref,alt = bases.split('>')
            
            positions.append(int(pos))

            # Add position after for indels
            if len(ref) != len(alt):
                positions.append(int(pos)+1)
                
            variant_data[sid]["variants"][var] = {}
            if var in sample["var_set"]:
                variant_data[sid]["variants"][var]["position"] = int(pos)
                if var in report_variants:
                    variant_data[sid]["variants"][var]["status"] = "match"
                else: 
                    variant_data[sid]["variants"][var]["status"] = "additional"
            else:
                variant_data[sid]["variants"][var]["status"] = "missing"

                
    depth_data_cursor = app.config['DEPTH_COLL'].find({'sample_oid': {'$in': sample_ids}, 'pos': {'$in': list(positions)}})
    depth_data = defaultdict(dict)
    for d in depth_data_cursor:
        depth_data[d["sample_oid"]][d["pos"]] = d
    
            
                    
    variant_annotation_cursor = app.config['VARIANT_COLL'].find({'_id': {'$in': list(variants)}})
    variant_annotation = {}

    
    for var_info in variant_annotation_cursor:
        csq = var_info["csq"]
        if "MUTATION" in csq and csq["MUTATION"] in variants_of_significance:
            var_info["significant"] = True
        elif "Protein_position" in csq and "SYMBOL" in csq and csq["SYMBOL"]+":"+str(csq["Protein_position"]) in positions_of_significance:
            var_info["significant_position"] = True
        variant_annotation[var_info["_id"]] = var_info

    return render_template('report.html', sample=sample_to_report, samples=samples_to_show, variant_data=variant_data, variants=variants, variant_annotation=variant_annotation, depth=depth_data, verbosity=verbosity, sample_counts=sample_counts, tot_samples=tot_samples, replicates=replicates, max_diff=max_diff)

@app.route('/mutations')
@login_required
def mutations():
    variant_annotations = app.config['VARIANT_COLL'].find()
    query = build_sample_query()
    samples = app.config['SAMPLE_COLL'].find(query)

    summarized_variants = defaultdict(dict)
    # summarize variants
    for sample in samples:
        if not "variants" in sample:
            continue
        
        for var in sample["variants"]:
            pos, change = var["id"].split("_")
            gene, aachange = var["aa"].split(":")

            if change in summarized_variants[pos]:
                summarized_variants[pos][change].append(sample)
            else:
                summarized_variants[pos][change] = [sample]

    return render_template('mutations.html', summary=summarized_variants, annotations=variant_annotations)
        
    
@app.route('/variant/<string:var_id>')
@login_required
def variant(var_id):
    variant_annotations = app.config['VARIANT_COLL'].find_one({'_id': unesc(var_id)})
    search_string = request.args.get("search_string", "")
    query = build_sample_query(search_string=search_string, mutation=unesc(var_id))

    samples_with_variant = app.config['SAMPLE_COLL'].find(query).sort('collection_date',-1)
    all_samples = app.config['SAMPLE_COLL'].find({'hidden':{'$ne':1}},{'collection_date':1, 'variants':1}).sort('collection_date')
    var_per_date = rolling_mean_variant(list(all_samples), var_id)
    samples_with_variant.rewind()
    
    return render_template('variant.html', annotations=variant_annotations, samples=samples_with_variant, var_per_date=var_per_date, var_id=var_id)


@app.route('/pangolin/<string:pango_type>')
@login_required
def pangolin(pango_type):
    search_string = request.args.get("search_string", "")
    query = build_sample_query(search_string=search_string, pango=pango_type)
    
    samples_with_type = app.config['SAMPLE_COLL'].find(query).sort('collection_date')
    samples = add_significant_variants(samples_with_type)

    all_samples = app.config['SAMPLE_COLL'].find({'hidden':{'$ne':1}},{'collection_date':1, 'pangolin.type':1}).sort('collection_date')
    pango_per_date, pango_total_counts = rolling_mean_pango_types(list(all_samples))
    
    lineages_of_concern = app.config["PANGO_LINEAGES_OF_CONCERN"]
    return render_template('pangolin.html', samples=samples, pangotype=pango_type, lineages_of_concern=lineages_of_concern, pango_per_date=pango_per_date)

@app.route('/bleek')
@login_required
def bleek():
    samples_with_type = app.config['SAMPLE_COLL'].find({'pangolin.type': "B.1.1.7", 'variants.id': '23012_G>A', 'hidden':{'$ne':1}}).sort('collection_date')
    samples = add_significant_variants(samples_with_type)

    lineages_of_concern = app.config["PANGO_LINEAGES_OF_CONCERN"]
    return render_template('bleek.html', samples=samples, lineages_of_concern=lineages_of_concern)


@app.route('/nextstrain/<string:clade>')
@login_required
def nextstrain(clade):
    clade = clade.replace("_", "/")
    samples_with_clade = app.config['SAMPLE_COLL'].find({'nextclade': clade, 'hidden':{'$ne':1}}).sort('collection_date')
    return render_template('nextstrain.html', samples=samples_with_clade, clade=clade)


@app.route('/createtree/<string:group>/<string:value>')
@login_required
def create_tree(group, value):

    if group == "pango":
        samples = app.config['SAMPLE_COLL'].find({'pangolin.type': value, 'hidden':{'$ne':1}, 'qc.pct_N_bases':{'$lte':app.config["QC_MAX_PCT_N"]}})
    if group == "nextstrain":
        value = value.replace('_','/')
        samples = app.config['SAMPLE_COLL'].find({'nextclade': value, 'hidden':{'$ne':1}, 'qc.pct_N_bases':{'$lte':app.config["QC_MAX_PCT_N"]}})
    if group == "bleek":
        samples = app.config['SAMPLE_COLL'].find({'pangolin.type': "B.1.1.7", 'variants.id':'23012_G>A', 'hidden':{'$ne':1}, 'qc.pct_N_bases':{'$lte':app.config["QC_MAX_PCT_N"]}})
    elif group == "all":
        samples = app.config['SAMPLE_COLL'].find({'hidden':{'$ne':1},'qc.pct_N_bases':{'$lte':app.config["QC_MAX_PCT_N"]}})

    presence = defaultdict(set)
    all_samples = set()
    sample_metadata = []
    pango_color = {}
    nextstrain_color = {}

    r = lambda: random.randint(0,255)
    
    for sample in samples:
        if "collection_date" not in sample:
            continue
        
        for var in sample.get("variants",[]):
            presence[var["id"]].add(sample["sample_id"])

        all_samples.add(sample["sample_id"])

        pango_type = sample["pangolin"]["type"]
        nextstrain_clade = sample.get("nextclade")

        # Generate random color for new pango types
        if pango_type not in pango_color:
            pango_color[pango_type] = '#%02X%02X%02X' % (r(),r(),r())

        if nextstrain_clade not in nextstrain_color:
            nextstrain_color[nextstrain_clade] = '#%02X%02X%02X' % (r(),r(),r())
        
        sample_metadata.append({
            'id': sample.get('sample_id'),
            'year': sample["collection_date"].year,
            'month': '{:02d}'.format(sample["collection_date"].month),
            'day': '{:02d}'.format(sample["collection_date"].day),
            'country': "Sweden",
            'pango': pango_type,
            'pango__color': pango_color[pango_type],
            'nextstrain': nextstrain_clade,
            'nextstrain__color': nextstrain_color[nextstrain_clade],
            'country__color': "#358"
        })
            
    distance = defaultdict(dict)

    for s1 in all_samples:
        for s2 in all_samples:
            distance[s1][s2] = 0
    
    for var_id, samples_with_variant in presence.items():
        samples_without_variant = all_samples ^ samples_with_variant
        for s_without in samples_without_variant:
            for s_with in samples_with_variant:
                distance[s_without][s_with] += 1
                distance[s_with][s_without] += 1

    data = []
    for s in all_samples:
        data.append(list(distance[s].values()))

    ids = list(all_samples)
    dm = DistanceMatrix(data, ids)
    tree = nj(dm)

    microreact_data = {
        "data": sample_metadata,
        "name": value,
        "tree": str(tree.root_at_midpoint()).rstrip(),
        "timeline_grouped": True,
    }

    response = requests.post('https://microreact.org/api/project/', headers={'Content-type': 'application/json; charset=UTF-8'}, data = json.dumps(microreact_data))

    if response.ok:
        return redirect(json.loads(response.content)["url"])

    return render_template('tree.html', tree=tree, response=response)
    


def get_similar_samples(sample_to_report, all_samples, max_diffs):

    # Collect variants in the sample to report
    report_variant_set = set()
    for var in sample_to_report.get("variants", []):
        report_variant_set.add(var["id"])

    sample_count_per_var = defaultdict(int)
        
    variants_of_similar_samples = report_variant_set.copy()

    sample_to_report["var_set"] = report_variant_set
    sample_to_report["n_diff"] = -1
    samples_to_show = [sample_to_report]
    samples_ids_to_show = [sample_to_report["_id"]]
    
    for sample_to_compare in all_samples:
        if sample_to_report["_id"] == sample_to_compare["_id"]:
            continue

        # Collect variants in the samples we're comparing to
        compare_variant_set = set()
        if "variants" in sample_to_compare:
            for var in sample_to_compare["variants"]:
                sample_count_per_var[var["id"]] += 1
                compare_variant_set.add(var["id"])

        # Calculate the number of differences in the variant sets
        n_diff = len(report_variant_set ^ compare_variant_set)

        # If samples are close save them
        if n_diff < max_diffs:
            variants_of_similar_samples |= compare_variant_set
            sample_to_compare["var_set"] = compare_variant_set
            sample_to_compare["n_diff"] = n_diff
            samples_to_show.append(sample_to_compare)

    return samples_to_show, variants_of_similar_samples, report_variant_set, sample_count_per_var


def add_significant_variants(samples):
    lineages_of_concern = app.config["PANGO_LINEAGES_OF_CONCERN"]
    variants_of_significance = set(app.config["VARIANTS_OF_BIOLOGICAL_SIGNIFICANCE"])
    positions_of_significance = set(app.config["POSITIONS_OF_BIOLOGICAL_SIGNIFICANCE"])
    
    samples_mod = []
    # Add significant variants
    for sample in samples:
        sample["significant_variants"] = []
        for variant in sample.get("variants", []):
            if variant["aa"] in variants_of_significance:
                sample["significant_variants"].append(variant)
        samples_mod.append(sample)

    return samples_mod


def criteria_is_selected(sample, exclude_criteria):
    sample_crit = sample.get("selection_criterion")
    if not exclude_criteria:
        return True
    if not exclude_criteria.get("surv") and sample_crit == "Allm&auml;n &ouml;vervakning":
        return False
    if not exclude_criteria.get("vaccine") and sample_crit == "Vaccingenombrott":
        return False
    if not exclude_criteria.get("reinfect") and sample_crit == "Reinfektion":
        return False
    if not exclude_criteria.get("travel") and sample_crit == "Reseanamnes":
        return False
    if not exclude_criteria.get("other") and sample_crit not in ["Vaccingenombrott", "Allm&auml;n &ouml;vervakning", "Reinfektion", "Reseanamnes"]:
        return False

    return True

def pango_per_week(samples, exclude_criteria=None):
    pango_count_per_week = defaultdict(lambda: defaultdict(int))
    weeks = defaultdict(int)
    
    for s in samples:
        if "collection_date" not in s or "pangolin" not in s or s["pangolin"]["type"] == "None":
            continue

        # Skip filtered selection criteria
        if not criteria_is_selected(s, exclude_criteria):
            continue
        
        pango = s["pangolin"]["type"]
        week_number = s["collection_date"].isocalendar()[1]
        weeks[week_number] += 1
        pango_count_per_week[pango]["Total"] += 1
        pango_count_per_week[pango][week_number] += 1

    return pango_count_per_week, weeks


def rolling_mean_pango_types(samples, exclude_criteria=None):
    pangolin_per_date = defaultdict(lambda: defaultdict(int))
    type_counts = defaultdict(int)
    for s in samples:
        if "collection_date" not in s or "pangolin" not in s or s["pangolin"]["type"] == "None":
            continue

        # Skip filtered selection criteria
        if not criteria_is_selected(s, exclude_criteria):
            continue        

        pango = s["pangolin"]["type"]
        type_counts[pango] += 1
        for i in range(7):
            day = s["collection_date"] + timedelta(i-3)
            pangolin_per_date[day.date()][pango] += 1

    return pangolin_per_date, type_counts


def rolling_mean_variant(samples, variant_id):
    variant_per_date = defaultdict(lambda: defaultdict(int))
    
    for s in samples:
        if "collection_date" not in s:
            continue

        for i in range(7):
            day = s["collection_date"] + timedelta(i-3)
            if "variants" in s:
                if any(v['id'] == variant_id for v in s["variants"]):
                    variant_per_date[day.date()][variant_id] += 1
                else:
                    variant_per_date[day.date()]['without_variant'] += 1

    return variant_per_date

def build_sample_query(search_string=None, hidden=False, validation=False, pango=None, mutation=None):
    query = {}

    if search_string:
        parts = search_string.split()
        query["$and"] = []
        for part in parts:
            query["$and"].append(
                {'$or': [
                    {'sample_id': {'$regex': part}},
                    {'pangolin.type': {'$regex':part}},
                    {'mlu': {'$regex': part}},
                    {'selection_criterion': {'$regex': part}},
                    {'variants.aa': {'$regex': part}},
                    {'lab': {'$regex': part}},
                    {'seqfacility': {'$regex': part}},
                ]})

    if pango:
        query["pangolin.type"] = pango

    if mutation:
        query["variants.id"] = mutation
            
    if not hidden:
        query["hidden"] = {'$ne':1}

    if not validation:
        query["validation"] = {'$ne':1}
            
    return query


@app.template_filter()
def one_letter_p(st):
    aa = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
          'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
          'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
          'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M', 'Ter': '*' }

    pattern = re.compile('|'.join(aa.keys()))
    new = pattern.sub(lambda x: aa[x.group()], st)
    return new

@app.template_filter()
def no_transid(nom):
    a = nom.split(':')
    if 1 < len(a):
        if a[1].startswith("c.") or a[1].startswith("p."):
            return a[1][2:]
        return a[1]
    
    return nom

@app.template_filter()
def unesc(st):
    if( len(st) > 0 ):
        return urllib.parse.unquote(st)
    else:
        return ""
                

@app.template_filter()
def pos_sort(l):
    return sorted(l, key=lambda k: int(k.split('_')[0]))

@app.template_filter()
def dictsort(l, key):
    return sorted(l, key=lambda k: k[key]) 

@app.template_filter()
def human_date(d):
    days_old = (datetime.now()-d).days

    # FIXME: This is weird and hardcoded. I can't find another way to get CET timestamps...
    cet = pytz.timezone("EET")
    d = d.astimezone(cet)
    
    if d.date() == datetime.today().date():
        return "Today "+d.strftime('%H:%M')
    yesterday = datetime.today().date() - timedelta(days=1)
    if d.date() == yesterday:
        return "Yesterday "+d.strftime('%H:%M')
    return d.strftime("%Y-%m-%d %H:%M")

@app.template_filter()
def date_notime(d):
    if d:
        return d.strftime("%Y-%m-%d")
    else :
        return ""


@app.template_filter()
def variant_link(v):
    covariants = {
        "S:A222V":"20A.EU1", "S:S477N":"20A.EU2", "S:N501Y":"S.N501", "S:E484K":"S.E484",
        "S:H69_V70del":"S.H69-", "S:N439K":"S.N439K", "S:Y453F":"S.Y453F", "S:S98F":"S.S98F",
        "S:L452R":"S.L452R", "S:D80Y":"S.D80Y", "S:A626S":"S.A626S", "S:V1122L":"S.V1122L",
        "S:Q677H":"S.Q677"
    }

    var_show = v
    if v == "S:H69_V70del":
        var_show = "S:H69-"
        
        
    if v in covariants:
        return f"<a href='https://covariants.org/variants/{ covariants[v] }'>{ var_show }</a>"

    return v

@app.template_filter()    
def pretty(value):
    return json.dumps(json.loads(value), sort_keys=True, indent=4, separators=(',', ': '))

@app.template_filter()
def get_dates(data):
    dates = []
    print(data)
    for date in sorted(data):
        tot = 0
        for i in data[date].values():
           tot += i
        if tot > 5:
            dates.append(date.strftime("%Y-%m-%d"))

    return dates

@app.template_filter()
def pct_type(data, key):
    pct = []
    for date in sorted(data):
        tot = 0
        for i in data[date].values():
           tot += i
        if tot > 5:
            selected_count =data[date].get(key, 0)
            pct.append(selected_count/tot)

    return pct

@app.template_filter()
def num_samples_per_date(data):
    num = []
    for date in sorted(data):
        tot = 0
        for i in data[date].values():
           tot += i
        if tot > 5:
            num.append(int(tot/7))

    return num


@login_manager.user_loader
def load_user(username):
    u = app.config['USERS_COLL'].find_one({"_id": username})
    if not u:
        return None
    return User(u['_id'], u['groups'], u.get('sc2_role', 'normal'))


@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if request.method == 'POST' and form.validate_on_submit():
        user = app.config['USERS_COLL'].find_one({"_id": form.username.data})
        if user and User.validate_login(user['password'], form.password.data):
            user_obj = User(user['_id'], user['groups'], user.get('sc2_role','normal'))
            login_user(user_obj)

            return redirect(request.args.get("next") or url_for("index"))

    return render_template('login.html', title='login', form=form)


@app.route('/logout')
def logout():
    logout_user()
    return redirect(url_for('login'))

