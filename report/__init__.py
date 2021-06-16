from flask import Flask
app = Flask(__name__)

import report.views

app.config.from_object('config')

