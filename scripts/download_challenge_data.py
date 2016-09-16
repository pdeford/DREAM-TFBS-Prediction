################################################################################
## Download competition data for the
## ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge
################################################################################

import shutil
import synapseclient
from synapseclient import Project, Folder, File
import sys

syn = synapseclient.Synapse()

# If you haven't set up a .synapseConfig file, you'll have to supply credentials
# syn.login(email = 'me@example.com', password = 'secret')
email = raw_input("Please enter Synapse Email: ")
password = raw_input("Password: ")
syn.login(email = email, password = password)

print "Make sure you've accepted the terms of use before running this script!"

# You may wish to copy these files to a specific destination directory. If so,
# set the path to that directory here or pass it as an argument to the script.
dest_dir = sys.argv[1] if len(sys.argv) > 1 else None

# -------------------------------------------------------------------------------
# All Challenge Data available for download:
# (* indicates data from the 'Essential Data Collection' - see https://www.synapse.org/#!Synapse:syn6131484/wiki/402033 )
# -------------------------------------------------------------------------------
# * ChIPseq fold_change_signal = syn6181334
# * ChIPseq labels = syn6181335
# * ChIPseq peaks conservative = syn6181337
# * ChIPseq peaks relaxed = syn6181338
# DNASE bams = syn6176232
# * DNASE fold_coverage_wiggles = syn6176233
# * DNASE peaks conservative = syn6176235
# * DNASE peaks relaxed = syn6176236
# * RNAseq = syn6176231
# * annotations = 'syn6184307'
# -------------------------------------------------------------------------------

# As written, this script will download the entire Essential Data Collection
# MODIFY THIS LINE TO INCLUDE THE SYNAPSE IDS OF DATA TYPES YOU WANT TO DOWNLOAD
folder_ids = ['syn6181334', 'syn6181335', 'syn6181337', 'syn6181338', 'syn6176233',
              'syn6176235', 'syn6176236', 'syn6176231', 'syn6184307']

data_files = []

for folder_id in folder_ids:

    # Get folder
    folder = syn.get(folder_id, downloadLocation=sys.argv[1])
    print 'Downloading contents of %s folder (%s)\n' % (folder.name, folder.id,)

    # Query for child entities
    query_results = syn.query('select id,name from file where parentId=="%s"' % folder_id)

    # Download all data files
    for entity in query_results['results']:
        print '\tDownloading file: ', entity['file.name']
        data_files.append(syn.get(entity['file.id'], downloadLocation=sys.argv[1]))

print 'Download complete!'


# Copy data files to dest_dir
"""
if dest_dir:
    print "Copying files to %s\n" % dest_dir
    for data_file in data_files:
        sys.stdout.write('.')
        shutil.copy2(data_file.path, dest_dir)

    print 'Copying complete!'
"""
