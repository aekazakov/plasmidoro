"""
    Annotate features with pLannotate
"""
import os
import uuid
import csv
import shutil
from pathlib import Path
from subprocess import Popen, PIPE, CalledProcessError
from amdplasmids.settings import TEMP_DIR
from magicpool.models import Feature, Feature_type

def export_plasmid_fasta(plasmid, working_dir):
    """
        Exports plasmid sequence in FASTA format. 
        Input:
            plasmid(Plasmid): Plasmid object
            working_dir(str): working directory, full path
        Output:
            path of the output FASTA file
    """
    out_file = os.path.join(working_dir, str(plasmid.id) + '.fasta')
    with open(out_file, 'w') as outfile:
        outfile.write('>' + str(plasmid.name) + '\n')
        outfile.write(str(plasmid.sequence))
    return out_file
    


def preprocess(plasmids, working_dir):
    """
        Creates all directories and input files. 
        Input:
            plasmids(list(Plasmid)): list of Plasmid objects
            working_dir(str): working directory, full path
        Output:
            path of the shell script running plannotate in conda environment
    """
    # Create directory
    if os.path.exists(working_dir) and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    
    # Create shell script
    plannotate_script = os.path.join(working_dir, 'run_plannotate.sh')
    with open(plannotate_script, 'w') as outfile:
        outfile.write('#!/usr/bin/bash\n')
        outfile.write('source /home/aekazakov/tools/anaconda3/etc/profile.d/conda.sh\n')
        outfile.write('conda activate plannotate\n')
        for plasmid in plasmids:
            if plasmid.sequence is None or plasmid.sequence == '':
                continue
                print(f"Plasmid {plasmid.name} has no sequence")
            outfile.write('echo "' + plasmid.name + '"\n')
            outfile.write(' '.join(['plannotate', 'batch', '-i',
                export_plasmid_fasta(plasmid, working_dir), '-c',
                '-o', working_dir
                ])
                + '\n')
        outfile.write('conda deactivate\n')
    return plannotate_script


def run(script_path):
    """
        Runs pLannotate for all plasmids.
    """
    cmd = ['/bin/bash', script_path]
    print(' '.join(cmd))
    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        # Suppress false positive no-member error (see https://github.com/PyCQA/pylint/issues/1860)
        # pylint: disable=no-member
        raise CalledProcessError(proc.returncode, proc.args)


def postprocess(plasmids, working_dir):
    """
        Finds pLannotate output files and creates new Features
        Input:
            plasmids(list(Plasmid)): list of Plasmid objects
            working_dir(str): working directory, full path
    """
    features_created = 0
    features_updated = 0
    for plasmid in plasmids:
        if plasmid.sequence is None or plasmid.sequence == '':
            continue
            print(f"Plasmid {plasmid.name} was skipped")
        print(plasmid.name)
        feature_type_dict = {item.name:item for item in Feature_type.objects.all()}
        in_file = os.path.join(working_dir, str(plasmid.id) + '_pLann.csv')
        with open(in_file, 'r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                feature_type_label = row['Type']
                feature_start = int(row['start location'])
                feature_end = int(row['end location'])
                feature_strand = int(row['strand'])
                feature_sequence = row['sequence']
                description = 'pLannotate: ' + row['Feature'] + ' [' + row['Type'] + '] - ' + row['percent identity'].split('.')[0] + '% ident.; ' + row['percent match length'].split('.')[0] + '% length; ' + row['Description']
                if feature_type_label not in feature_type_dict:
                    Feature_type.objects.create(name=feature_type_label)
                    feature_type_dict = {item.name:item for item in Feature_type.objects.all()}
                feature_type = feature_type_dict[feature_type_label]
                if Feature.objects.filter(plasmid = plasmid,
                        start = feature_start,
                        end = feature_end,
                        strand = feature_strand,
                        feature_type = feature_type
                        ).exists():
                    # Feature exists. Update description.
                    feature = Feature.objects.get(plasmid = plasmid,
                        start = feature_start,
                        end = feature_end,
                        strand = feature_strand,
                        feature_type = feature_type
                        )
                    if description not in feature.description:
                        feature.description = feature.description + ' ' + description
                        feature.save()
                        features_updated += 1
                    continue
                else:
                    # Feature does not exist. Create it.
                    if feature_strand == -1:
                        location_str = f'[{feature_start}:{feature_end}](-)'
                    else:
                        location_str = f'[{feature_start}:{feature_end}](+)'
                        feature = Feature.objects.create(name = row['Feature'].replace(' ', '_'),
                            plasmid = plasmid,
                            feature_type = feature_type,
                            sequence = row['sequence'],
                            sequence_id = plasmid.name,
                            start = feature_start,
                            end = feature_end,
                            strand = feature_strand,
                            location_str = location_str,
                            description = description
                            )
                        features_created += 1
    print('Features created: ', features_created)
    print('Features updated: ', features_updated)
    _cleanup(working_dir)


def _cleanup(working_dir):
    shutil.rmtree(working_dir)

def application(plasmids):
    """
        This function is an entry point of the module.
        Input:
            plasmids(list(Plasmid)): list of Plasmid objects
    """
    working_dir = os.path.join(TEMP_DIR,str(uuid.uuid4()))
    
    script_path = preprocess(plasmids, working_dir)
    run(script_path)
    postprocess(plasmids, working_dir)
    
