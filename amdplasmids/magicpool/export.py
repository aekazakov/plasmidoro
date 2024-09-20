"""
    Export to text or xlsx files
"""
import csv
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
    
