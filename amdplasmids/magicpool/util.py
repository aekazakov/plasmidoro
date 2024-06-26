"""
    Various utility functions
"""
import os
import sys
import gzip
import openpyxl
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio import GenBank
from magicpool.models import *

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def guess_feature_type(name, type_label):
    feature_type_dict = {item.name:item for item in Feature_type.objects.all()}

    if type_label in feature_type_dict:
        return feature_type_dict[type_label]
    elif type_label == 'CDS':
        return Feature_type.objects.get(name='gene')
    elif name.endswith('romoter'):
        return Feature_type.objects.get(name='promoter')
    elif 'origin of replication' in name:
        return Feature_type.objects.get(name='origin')
    elif name.startswith('Ori'):
        return Feature_type.objects.get(name='origin')
    elif name == 'ColE1':
        return Feature_type.objects.get(name='origin')
    elif 'BsmBI' in name:
        return Feature_type.objects.get(name='restriction site')
    elif 'inverted repeat' in name:
        return Feature_type.objects.get(name='IR')
    elif 'transposase enzyme' in name:
        return Feature_type.objects.get(name='gene')
    else:
        return Feature_type.objects.get(name='unknown')

def create_feature(feature, plasmid_obj):
    feature_type_label = feature.type
    if feature_type_label == 'source':
        return 0
    feature_name = ''
    if 'locus_tag' in feature.qualifiers:
        feature_name = str(feature.qualifiers['locus_tag'][0])
    elif 'label' in feature.qualifiers:
        feature_name = str(feature.qualifiers['label'][0])
    if feature_name == '':
        print('This feature has no name:')
        print(feature)
        return 0
    description = ''
    if 'gene' in feature.qualifiers:
        description += str(';' .join(feature.qualifiers['gene']))
    if 'label' in feature.qualifiers:
        if description != '':
            description += '; '
        description += str(';' .join(feature.qualifiers['label']))
    if 'note' in feature.qualifiers:
        if description != '':
            description += '; '
        description += str(feature.qualifiers['note'][0])
    feature_type = guess_feature_type(feature_name, feature_type_label)
    feature_start = int(feature.location.start)
    feature_end = int(feature.location.end)
    feature_strand = int(feature.location.strand)
    feature_obj = Feature.objects.create(
        name = feature_name,
        plasmid = plasmid_obj,
        feature_type = feature_type,
        sequence = plasmid_obj.sequence[feature_start:feature_end],
        sequence_id = plasmid_obj.name,
        start = feature_start,
        end = feature_end,
        strand = feature_strand,
        location_str = str(feature.location),
        description = description
        )
    if feature_type.name == 'gene' or 'translation' in feature.qualifiers:
        if 'gene' in feature.qualifiers:
            protein_name = str(feature.qualifiers['gene'][0])
        elif 'label' in feature.qualifiers:
            protein_name = str(feature.qualifiers['label'][0])
        elif 'locus_tag' in feature.qualifiers:
            protein_name = str(feature.qualifiers['locus_tag'][0])
        else:
            protein_name = 'unknown_protein'

        if 'product' in feature.qualifiers:
            protein_function = str(feature.qualifiers['product'][0])
        elif 'note' in feature.qualifiers:
            protein_function = str(feature.qualifiers['note'][0])
        elif 'label' in feature.qualifiers:
            protein_function = str(feature.qualifiers['label'][0])
        else:
            protein_function = 'unidentified function'

        if 'translation' in feature.qualifiers:
            protein_sequence = ''.join(feature.qualifiers['translation'])
        else:
            protein_sequence = ''
            
        protein_obj = Protein.objects.create(
            name = protein_name,
            sequence = protein_sequence,
            function = protein_function,
            feature = feature_obj
        )
    return 1    
            

def import_seq_records(records, seq_name, overwrite_existing=False):
    record_count = 0
    for seq_record in records:
        if record_count > 1:
            continue
        existing_plasmids = Plasmid.objects.filter(name=seq_name)
        seq_record_sequence = str(seq_record.seq)
        feature_count = 0
        if existing_plasmids.exists():
            # Check sequence
            existing_plasmid = existing_plasmids.first()
            if not overwrite_existing:
                print(existing_plasmid.name, 'object not updated because the --overwrite option is omitted')
                ret_val = 0
            elif existing_plasmid.sequence == seq_record_sequence and Feature.objects.filter(plasmid=existing_plasmid.id).exists:
                # Nothing to do here
                print(existing_plasmid.name, 'object not updated because the sequence did not change')
                ret_val = 0
            else:
                # Replace sequence, wipe out old features and create new features
                existing_plasmid.sequence = seq_record_sequence
                existing_plasmid.save()
                Feature.objects.filter(plasmid=existing_plasmid.id).delete()
                Protein.objects.filter(feature=None).delete()
                
                for feature in seq_record.features:
                    feature_count += create_feature(feature, existing_plasmid)
                print(existing_plasmid.name, 'object updated')
                ret_val = 1
        else:
            # Create new plasmid object and features
            plasmid_obj = Plasmid.objects.create(
                name = seq_name,
                amd_number = '',
                description = '',
                sequence = seq_record_sequence,
            )
            for feature in seq_record.features:
                feature_count += create_feature(feature, plasmid_obj)
            print(seq_name, ' object created')
            ret_val = 1
        record_count += 1
    if record_count > 1:
        print(record_count, 'records provided for import. Only the first record will be imported')
    print(feature_count, 'features created')
    return ret_val

def import_plasmid_snapgene(dna_file, overwrite_existing):
    print('Working on Snapgene file', dna_file)
    records = SeqIO.parse(dna_file, "snapgene")
    filename = dna_file.split('/')[-1]
    ret_val = import_seq_records(records, filename.split('.dna')[0], overwrite_existing)
    return ret_val
            
def import_plasmid_gbk(gbk_file, overwrite_existing):
    print('Working on GenBank file', gbk_file)
    if gbk_file.endswith('.gz'):
        gbk_handle = gzip.open(gbk_file, 'rt')
    else:
        gbk_handle = open(gbk_file, 'r')
    records = SeqIO.parse(gbk_handle, "genbank")
    filename = gbk_file.split('/')[-1]
    ret_val = import_seq_records(records, '.'.join(filename.split('.')[:-1]), overwrite_existing)
    gbk_handle.close()
    return ret_val

def update_plasmids(work_dir, overwrite_existing):
    obj_count = 0
    for root, subdirs, files in os.walk(work_dir):
        for filename in files:
            if filename.endswith('.dna'):
                obj_count += import_plasmid_snapgene(os.path.join(root, filename), overwrite_existing)
            if filename.endswith('.gb') or filename.endswith('.ape'):
                obj_count += import_plasmid_gbk(os.path.join(root, filename), overwrite_existing)
    print(obj_count, 'plasmids created and/or updated')

    
def create_plasmid(plasmid_name, plasmid_data):
    amd_number = ''
    description = ''
    if 'AMD number' in plasmid_data:
        amd_number = plasmid_data['AMD number']
    if 'Description' in plasmid_data:
        description = plasmid_data['Description']
    plasmid_obj = Plasmid.objects.create(
        name = plasmid_name,
        amd_number = amd_number,
        description = description
    )
        
    if 'Drug marker' in plasmid_data:
        for marker in plasmid_data['Drug marker'].split('; '):
            try:
                marker_obj = Drug_marker.objects.get(name=marker)
            except Drug_marker.DoesNotExist:
                marker_obj = Drug_marker.objects.create(
                    name = marker,
                    drug = marker,
                    note = ''
                    )
            plasmid_obj.drug_markers.add(marker_obj)
            plasmid_obj.save()
    for key, value in plasmid_data.items():
        if key in ('AMD number', 'Description', 'Drug marker'):
            continue
        plasmid_info_obj = Plasmid_info.objects.create(
            plasmid = plasmid_obj,
            param = key,
            value = value
        )
    return 1
    
def import_plasmids_table(xlsx_path):
    xlsx_path = Path(xlsx_path)
    wb_obj = openpyxl.load_workbook(xlsx_path)
    sheet = wb_obj.active
    xlsx_header = []
    data_imported = defaultdict(dict)
    created_count = 0
    for i, row in enumerate(sheet.iter_rows(values_only=True)):
        if i == 0:
            xlsx_header = row[1:]
        else:
            plasmid_name = row[0]
            if plasmid_name == '' or plasmid_name is None:
                continue
            for j, cell in enumerate(row[1:]):
                if cell != '' and cell != 'None' and cell is not None:
                    data_imported[plasmid_name][xlsx_header[j]] = str(cell)
    for plasmid_name, plasmid_data in data_imported.items():
        try:
            plasmid = Plasmid.objects.get(name=plasmid_name)
            for key,value in plasmid_data.items():
                if key == 'AMD number':
                    if plasmid.amd_number != value:
                        plasmid.amd_number = value
                        plasmid.save()
                elif key == 'Drug marker':
                    existing_markers = [item.name for item in plasmid.drug_markers.all()]
                    for marker in value.split('; '):
                        if marker not in existing_markers:
                            try:
                                marker_obj = Drug_marker.objects.get(name=marker)
                            except Drug_marker.DoesNotExist:
                                marker_obj = Drug_marker.objects.create(
                                    name = marker,
                                    drug = marker,
                                    note = ''
                                    )
                            plasmid.drug_markers.add(marker_obj)
                            plasmid.save()
                elif key == 'Description':
                    if value not in plasmid.description:
                        plasmid.description += value
                        plasmid.save()
                else:
                    if Plasmid_info.objects.filter(plasmid=plasmid.id,param=key).exists():
                        plasmid_info_obj = Plasmid_info.objects.get(plasmid=plasmid.id,param=key)
                        if plasmid_info_obj.value != value:
                            plasmid_info_obj.value = value
                            plasmid_info_obj.save()
                    else:
                        Plasmid_info.objects.create(
                            plasmid = plasmid,
                            param = key,
                            value = value
                            )
        except Plasmid.DoesNotExist:
            created_count += create_plasmid(plasmid_name, plasmid_data)
    print(created_count, 'new plasmids created')

def import_magic_pool(xlsx_path):
    xlsx_path = Path(xlsx_path)
    wb_obj = openpyxl.load_workbook(xlsx_path, data_only = True)
    sheet = wb_obj.get_sheet_by_name('part_names_overlaps')
    xlsx_header = []
    created_count = 0
    existing_parts = {item.name:item for item in Magic_pool_part.objects.all()}
    for i, row in enumerate(sheet.iter_rows()):
        if i == 0:
            xlsx_header = [cell.value for cell in row[1:]]
        else:
            part_name = row[0].value
            if part_name == '' or part_name is None:
                continue
            if part_name in existing_parts:
                continue
            existing_overhangs = {item.name:item for item in Overhang.objects.all()}
            description = ''
            upstream_overhang = None
            upstream_overhang_color = None
            downstream_overhang = None
            downstream_overhang_color = None
            partinfo = []
            for j, cell in enumerate(row[1:]):
                if xlsx_header[j] == '' or xlsx_header[j] is None:
                    break
                cell_text = cell.value
                if cell_text != '' and cell_text != 'None' and cell_text is not None:
                    if xlsx_header[j] == 'Contains':
                        description = cell_text
                    elif xlsx_header[j] == '4bp overhang upstream':
                        upstream_overhang = cell_text
                        cell_color = int(cell.fill.start_color.index, 16)
                        upstream_overhang_color = "%06x" % (cell_color and 0xFFFFFF)
                    elif xlsx_header[j] == '4bp overhang downstream':
                        downstream_overhang = cell_text
                        cell_color = int(cell.fill.start_color.index, 16)
                        downstream_overhang_color = "%06x" % (cell_color and 0xFFFFFF)
                    else:
                        partinfo.append((xlsx_header[j], cell_text))
            if upstream_overhang is None:
                print(part_name, ': Upstream overhang not found')
                continue
            if downstream_overhang is None:
                print(part_name, ': Downstream overhang not found')
                continue
            if upstream_overhang in existing_overhangs:
                upstream_overhang_obj = existing_overhangs[upstream_overhang]
            else:
                upstream_overhang_obj = Overhang.objects.create(
                    name = upstream_overhang,
                    sequence = upstream_overhang,
                    color = upstream_overhang_color
                    )
            if downstream_overhang in existing_overhangs:
                downstream_overhang_obj = existing_overhangs[downstream_overhang]
            else:
                downstream_overhang_obj = Overhang.objects.create(
                    name = downstream_overhang,
                    sequence = downstream_overhang,
                    color = downstream_overhang_color
                    )
            part_obj = Magic_pool_part.objects.create(
                name = part_name,
                description = description,
                upstream_overhang = upstream_overhang_obj,
                downstream_overhang = downstream_overhang_obj,
                )
            for item in partinfo:
                info_obj = Magic_pool_part_info.objects.create(
                    magic_pool_part = part_obj,
                    name = item[0],
                    description = item[1]
                )
            created_count += 1

    sheet = wb_obj.get_sheet_by_name('overview')
    xlsx_header = []
    existing_vectors = {item.name:item for item in Vector.objects.all()}
    current_vector_name = None
    current_description = ''
    magic_pool_parts = []
    for i, row in enumerate(sheet.iter_rows(values_only=True)):
        if i == 0:
            continue
        else:
            cell_a = str(row[0])
            cell_b = str(row[1])
            cell_c = str(row[2])
            if cell_a != '' and cell_a != 'None' and cell_a is not None:
                if current_vector_name is not None and current_vector_name not in existing_vectors:
                        vector_obj = Vector.objects.create(
                            name = current_vector_name,
                            description = current_description
                        )
                        for item_index, item in enumerate(magic_pool_parts):
                            try:
                                magic_pool_part_obj = Magic_pool_part.objects.get(name=item)
                                vector_part_obj = Vector_part.objects.create(
                                    vector = vector_obj,
                                    part = magic_pool_part_obj,
                                    order = item_index
                                )
                            except Magic_pool_part.DoesNotExist:
                                print(item, ' NOT FOUND IN MAGIC POOL PARTS')
                magic_pool_parts = []
                current_vector_name = cell_a
                current_description = cell_b
                magic_pool_parts.append(cell_c)
            elif cell_c != '' and cell_c != 'None' and cell_c is not None:
                magic_pool_parts.append(cell_c)
    if magic_pool_parts and current_vector_name is not None:
        if current_vector_name not in existing_vectors:
            vector_obj = Vector.objects.create(
                name = current_vector_name,
                description = current_description
            )
            for item_index, item in enumerate(magic_pool_parts):
                try:
                    magic_pool_part_obj = Magic_pool_part.objects.get(name=item)
                    vector_part_obj = Vector_part.objects.create(
                        vector = vector_obj,
                        part = magic_pool_part_obj,
                        order = item_index
                    )
                except Magic_pool_part.DoesNotExist:
                    print(item, ' NOT FOUND IN MAGIC POOL PARTS')
                
                    
    
