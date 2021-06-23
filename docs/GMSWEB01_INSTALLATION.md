# SC2REPORTER CONFIG

This is more or less just a log of what I did to install sc2reporter on the GMSWEB01 virtual machine.

# Preparatory setup

## Install some basic requirements
```
# dnf groupinstall "Development Tools"
# dnf install git python38-mod_wsgi python38-devel perl-App-cpanminus
# cpanm MongoDB
```

## Install and enable apache
```
# dnf install httpd
# systemctl enable httpd
# systemctl start httpd
```

## Open firewall
```
# firewall-cmd --zone=public --permanent --add-service=http
# firewall-cmd --reload
```

## Install and configure mongodb 

create `/etc/yum.repos.d/mongodb-org-4.4.repo`, with this content:

```
[mongodb-org-4.4]
name=MongoDB Repository
baseurl=https://repo.mongodb.org/yum/redhat/$releasever/mongodb-org/4.4/x86_64/
gpgcheck=1
enabled=1
gpgkey=https://www.mongodb.org/static/pgp/server-4.4.asc
```

Then run:
```
# dnf install -y mongodb-org
```

### Make sure mongodb plays nice with SELinux
```
$ cat > mongodb_cgroup_memory.te <<EOF
module mongodb_cgroup_memory 1.0;
require {
      type cgroup_t;
      type mongod_t;
      class dir search;
      class file { getattr open read };
}
#============= mongod_t ==============
allow mongod_t cgroup_t:dir search;
allow mongod_t cgroup_t:file { getattr open read };
EOF
$ checkmodule -M -m -o mongodb_cgroup_memory.mod mongodb_cgroup_memory.te
$ semodule_package -o mongodb_cgroup_memory.pp -m mongodb_cgroup_memory.mod
$ sudo semodule -i mongodb_cgroup_memory.pp
$ sudo /usr/sbin/setsebool -P httpd_can_network_connect 1
```

### Start mongod
```
$ sudo systemctl start mongod
```



## Install sc2reporter

### Create directory and give it reasonable permissions
```
# mkdir /data/sc2reporter
# chown xhabjo:adm-gms-submit01-p-ssh@GU.GU.SE /data/sc2reporter/
# cd /data/sc2reporter
```

### Clone repo
```
$ git clone https://github.com/bjhall/sc2reporter.git
$ cd sc2reporter
```

### Set up virtual environment and install requirements
```
$ python3.8 -m venv venv
$ source venv/bin/activate 
$ pip install --upgrade pip
$ pip install -r requirements.txt
```

### Add apache configuration and restart apache
$ sudo cp sc2reporter.conf /etc/httpd/conf.d/
$ sudo apachectl restart
