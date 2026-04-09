"""
    Export to text or xlsx files
"""
import csv
import gzip
from io import BytesIO
from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from django.http import HttpResponse
from magicpool.models import *


def export_plasmids(request):
    # Create the HttpResponse object with the appropriate CSV header.
    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="exported_plasmids.tab"'
    writer = csv.writer(response, delimiter='\t')
    query = request.GET.get('query')
    query_type = request.GET.get('type')
    if query_type == 'vector_id':
        vector = Vector_type.objects.get(id=query)
        magic_pool_part_type_ids = Magic_pool_part_type.objects.filter(vector_type_part__vector_type__id=query).values_list('id')
        plasmids = Plasmid.objects.filter(magic_pool_part__in=magic_pool_part_type_ids)
    else:
        plasmids = Plasmid.objects.none()
    writer.writerow([
                     'Name',
                     'AMD number',
                     'Designation',
                     'Magic pool part',
                     'Description',
                     ])
    for item in plasmids:
        writer.writerow([item.name,
                         item.amd_number,
                         item.magic_pool_designation,
                         item.magic_pool_part.name,
                         item.description
                         ])
    return response


def export_plasmid(plasmid, output_buffer):
    '''
        Reads one genome from GenBank file,
        adds gene mappings and annotations from the GenomeDepot database,
        writes the genome in GenBank format into the output buffer
    '''
    seqrecord = SeqRecord(
        Seq(str(plasmid.sequence)),
        id=plasmid.amd_number,
        name=plasmid.name,
        annotations = {'molecule_type' : 'DNA', 'accessions':[]}
    )
    seqrecord.description = 'Plasmid ' + plasmid.name + ', complete sequence'
    seqrecord.annotations['comment'] = plasmid.description
    seqrecord.annotations['organism'] = 'Plasmid ' + plasmid.name
    seqrecord.annotations['source'] = 'Plasmid ' + plasmid.name
    seqrecord.annotations['accessions'].append(plasmid.amd_number)
    seqrecord.annotations['references'] = []
    ref = SeqFeature.Reference()
    ref.title = 'Deutschbauer Lab plasmids'
    ref.comment = 'https://iseq.lbl.gov/plasmids/plasmid/' + str(plasmid.id)
    seqrecord.annotations['references'].append(ref)
    if plasmid.sequence_file:
        ref = SeqFeature.Reference()
        ref.title = 'Deutschbauer Lab Google Drive'
        ref.comment = str(plasmid.sequence_file)
        seqrecord.annotations['references'].append(ref)
    source_location = sf.FeatureLocation(sf.ExactPosition(0), sf.ExactPosition(len(seqrecord)), strand=+1)
    source_feature = sf.SeqFeature(source_location, type='source', qualifiers={'organism':plasmid.name,'mol_type':'genomic DNA', 'db_xref':'taxon:29278'})
    seqrecord.features.append(source_feature)
    for feature in Feature.objects.filter(plasmid=plasmid):
        feature_location = sf.FeatureLocation(feature.start, feature.end, strand=feature.strand)
        feature_obj = sf.SeqFeature(feature_location, type=feature.feature_type.name, qualifiers={'note':feature.description,'name':feature.name,'location':feature.location_str})
        seqrecord.features.append(feature_obj)
    seqrecord.features = sorted(seqrecord.features,
                                key=lambda x: x.location.start
                                )
    SeqIO.write(seqrecord, output_buffer, 'genbank')    


def export_gbk(request, name):
    response = HttpResponse(content_type='application/x-gzip')
    response['Content-Disposition'] = f'attachment; filename={name}.gbk'
    response['Content-Encoding'] = 'gzip'
    try:
        plasmid = Plasmid.objects.get(name = name)
    except Genome.DoesNotExist:
        logger.error('Genome not found: ' + str(name))
        return render(request,
              '404.html',
              {'searchcontext': 'Plasmid ' + name + ' does not exist'}
              )
    gzip_buffer = BytesIO()
    with gzip.open(gzip_buffer, 'wt') as outfile:
        export_plasmid(plasmid, outfile)
    response.write(gzip_buffer.getvalue())
    return response
