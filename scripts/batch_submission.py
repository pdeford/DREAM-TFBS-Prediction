#!/usr/bin/env python 
# 
# Submit files to challenge
# 
doc = """Usage: make-submissions SUBMISSIONFILES..."""

from docopt import docopt
import synapseclient
from synapseclient import File

# 
# Parse options
# 
opts = docopt(doc) 
submission_files = opts['SUBMISSIONFILES']

# 
# Log in 
# 
print('Logging in.')
syn = synapseclient.Synapse()
email = raw_input("Please enter Synapse Email: ")
password = raw_input("Password: ")
syn.login(email = email, password = password)

# 
# Our project 
# 
print('Getting project.')
project = syn.get(7114158) 

# 
# Get the object representing the evaluation queue
# 
print('Getting evaluation queue.')
# Ladder Round 
evaluation = syn.getEvaluation(7071644)
# Final Round
# evaluation = syn.getEvaluation(7212779)


for filename in submission_files: 
    #
    # Upload our file 
    #
    print('Uploading: {}'.format(filename)) 
    upload_file = File(filename, description='Submission', parentId = project.id) 
    upload_file = syn.store(upload_file)

    #
    # Submit the file to the evaluation queue
    #
    print('Submitting: {}'.format(filename))
    submission = syn.submit(evaluation, upload_file, name='Submission', team='ChIP Shape')