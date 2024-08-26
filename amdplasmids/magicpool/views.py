import os
from django.shortcuts import render
from django.template import loader
from django.http import HttpResponse
from django.db.models import Q
from magicpool.models import *
from magicpool.blast_search import run_protein_search
from magicpool.blast_search import run_nucleotide_search
from amdplasmids.settings import STATICFILES_DIRS, STATIC_URL

def index(request):
	#return HttpResponse("This is a start page.")
    template = loader.get_template('magicpool/index.html')
    context = {'site_title':'Plasmidoro'}
    return HttpResponse(template.render(context, request))


def show_help(request):
    '''
        Displays help page
    '''
    template = loader.get_template('magicpool/help.html')
    context = {}
    return HttpResponse(template.render(context, request))


def plasmid_detail(request, plasmid_id):
    '''
        Displays plasmid page
    '''
    template = loader.get_template('magicpool/plasmid.html')
    plasmid = Plasmid.objects.get(id=plasmid_id)
    context = {
        'site_title':plasmid.name + ' [' + plasmid.amd_number + ']',
        'plasmid':plasmid,
        'plasmid_info':Plasmid_info.objects.filter(plasmid=plasmid.id),
        'features': Feature.objects.filter(plasmid=plasmid.id).order_by('start')
        
        }
    if Strain.objects.filter(amd_number=plasmid.amd_number).exists():
        context['strain'] = Strain.objects.get(amd_number=plasmid.amd_number)
    print(context)
    return HttpResponse(template.render(context, request))


def plasmid_viewer(request, plasmid_id):
    '''
        Displays plasmid viewer
    '''
    template = loader.get_template('magicpool/viewer.html')
    plasmid = Plasmid.objects.get(id=plasmid_id)
    
    context = {
        'site_title':plasmid.name + ' [' + plasmid.amd_number + '] viewer',
        'plasmid':plasmid,
        'features':Feature.objects.filter(plasmid=plasmid_id),
        'translations':Feature.objects.filter(plasmid=plasmid_id, feature_type__name='gene'),
        }
    print(context)
    return HttpResponse(template.render(context, request))


def part_detail(request, part_id):
    '''
        Displays magic pool part page
    '''
    template = loader.get_template('magicpool/part.html')
    part = Magic_pool_part.objects.get(id=part_id)
    context = {
        'site_title':part.name,
        'part':part,
        'info':Magic_pool_part_info.objects.filter(magic_pool_part=part_id),
        'vectors':Vector.objects.filter(vector_part__part__id=part_id),
        'plasmids':Plasmid.objects.filter(magic_pool_part__id=part_id)
        }
    return HttpResponse(template.render(context, request))


def vector_detail(request, vector_id):
    '''
        Displays vector design page
    '''
    template = loader.get_template('magicpool/vector.html')
    vector = Vector.objects.get(id=vector_id)
    img_path = os.path.join(STATICFILES_DIRS[0], 'vectors', vector.name + '.png')
    context = {'vector':vector,
        'info':Vector_info.objects.filter(vector=vector_id),
        'vector_parts':Vector_part.objects.filter(vector=vector_id).select_related('part'),
        }
    if os.path.exists(img_path):
        context['img'] = STATIC_URL + '/vectors/' + vector.name + '.png'
    return HttpResponse(template.render(context, request))

    
def oligo_detail(request, oligo_id):
    '''
        Displays oligo page
    '''
    template = loader.get_template('magicpool/oligo.html')
    oligo = Oligo.objects.get(id=oligo_id)
    context = {'oligo':oligo,
        'info':Oligo_info.objects.filter(oligo=oligo_id),
        }
    return HttpResponse(template.render(context, request))
    

def strain_detail(request, strain_id):
    '''
        Displays strain page
    '''
    template = loader.get_template('magicpool/strain.html')
    strain = Strain.objects.get(id=strain_id)
    plasmids = Plasmid.objects.filter(name=strain.plasmid)
    context = {
        'site_title': str(strain),
        'strain':strain,
        'strain_info':Strain_info.objects.filter(strain=strain.id),
        'plasmids':plasmids
        }
        
    print(context)
    return HttpResponse(template.render(context, request))

    
def plasmids(request):
    '''
        Displays list of plasmids
    '''
    template = loader.get_template('magicpool/plasmids.html')
    context = {
        'itemlist':Plasmid.objects.order_by('name')
        }
    return HttpResponse(template.render(context, request))


def parts(request):
    '''
        Displays list of parts
    '''
    template = loader.get_template('magicpool/parts.html')
    context = {
        'itemlist':Magic_pool_part.objects.order_by('name')
        }
    return HttpResponse(template.render(context, request))


def vectors(request):
    '''
        Displays list of vectors
    '''
    template = loader.get_template('magicpool/vectors.html')
    context = {
        'itemlist':Vector.objects.order_by('name')
        }
    print(context)
    return HttpResponse(template.render(context, request))


def strains(request):
    '''
        Displays list of strains
    '''
    template = loader.get_template('magicpool/strains.html')
    context = {
        'itemlist':Strain.objects.order_by('amd_number')
        }
    print(context)
    return HttpResponse(template.render(context, request))


def nucleotidesearchform(request):
    '''
        Displays nucleotide sequence search form
    '''
    return render(request,'magicpool/nucleotidesearchform.html')


def nucleotidesearch(request):
    '''
        Returns nucleotide search results 
    '''
    context = {}
    if request.POST.get("sequence"):
        result = {}
        params = {'sequence': request.POST.get("sequence"),
                  'evalue': request.POST.get("evalue"),
                  'hitstoshow': request.POST.get("hitstoshow")
                  }
        hits, searchcontext, query_len, _ = run_nucleotide_search(params)
        if searchcontext != '':
            context['searchcontext'] = searchcontext
        for row in hits:
            row=row.split('\t')
            print('Search for gene %s', row[1])
            if row[0] not in result:
                result[row[0]] = []
            unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
            query_cov = (query_len - unaligned_part) * 100.0 / query_len
            target_tokens = row[1].split('|')
            target_id = int(target_tokens[0])
            if target_tokens[2] == 'plasmid':
                plasmid = Plasmid.objects.get(id = target_id)
                if plasmid.amd_number != '':
                    plasmid_label = plasmid.name + '/' + plasmid.amd_number
                else:
                    plasmid_label = plasmid.name
                hit = ['Plasmid',
                       plasmid.id,
                       plasmid_label,
                       '{:.1f}'.format(float(row[2])),
                       row[3],
                       '{:.1f}'.format(query_cov),
                       row[10],
                       row[11]
                       ]
                result[row[0]].append(hit)
            elif target_tokens[2] == 'oligo':
                oligo = Oligo.objects.get(id = target_id)
                hit = ['Oligo',
                       oligo.id,
                       oligo.name,
                       '{:.1f}'.format(float(row[2])),
                       row[3],
                       '{:.1f}'.format(query_cov),
                       row[10],
                       row[11]
                       ]
                result[row[0]].append(hit)
        context['searchresult'] = result
    else:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Provide nucleotide sequence in FASTA format'}
                      )
    return render(request, 'magicpool/nucleotidesearch.html', context)

    '''
        Displays nucleotide sequence search form
    '''
    #return render(request,'magicpool/nucleotidesearchform.html')
    return HttpResponse("Not implemented.")


def proteinsearchform(request):
    '''
        Displays protein sequence search form
    '''
    return render(request,'magicpool/proteinsearchform.html')


def proteinsearch(request):
    '''
        Returns protein search results 
    '''
    context = {}
    if request.POST.get("sequence"):
        result = {}
        params = {'sequence': request.POST.get("sequence"),
                  'evalue': request.POST.get("evalue"),
                  'hitstoshow': request.POST.get("hitstoshow")
                  }
        hits, searchcontext, query_len, _ = run_protein_search(params)
        if searchcontext != '':
            context['searchcontext'] = searchcontext
        for row in hits:
            row=row.split('\t')
            print('Search for gene %s', row[1])
            if row[0] not in result:
                result[row[0]] = []
            unaligned_part =  int(row[6]) - 1 + query_len - int(row[7])
            query_cov = (query_len - unaligned_part) * 100.0 / query_len
            proteins = Protein.objects.select_related(
                'feature', 'feature__plasmid'
            ).filter(
                id = int(row[1].split('|')[0])
            )
            for protein in proteins:
                hit = [protein.name,
                       protein.feature.plasmid.id,
                       protein.feature.plasmid.name,
                       protein.function,
                       '{:.1f}'.format(float(row[2])),
                       row[3],
                       '{:.1f}'.format(query_cov),
                       row[10],
                       row[11]
                       ]
                result[row[0]].append(hit)
        context['searchresult'] = result
    else:
        return render(request,
                      '404.html',
                      {'searchcontext': 'Provide protein sequence in FASTA format'}
                      )
    return render(request, 'magicpool/proteinsearch.html', context)


def textsearchform(request):
    '''
        Displays text search form
    '''
    return render(request,'magicpool/search.html')


def textsearch(request):
    '''
        Search in plasmids
    '''
    query = request.GET.get('query')
    type = request.GET.get('type')
    if type == 'plasmid':
        if query:
            object_list = Plasmid.objects.filter(
                    (Q(name__icontains=query) |
                    Q(amd_number__icontains=query) |
                    Q(description__icontains=query))
                ).distinct().order_by('name')
        else:
            query = ''
            object_list = Plasmid.objects.none()
        template = loader.get_template('magicpool/plasmids.html')
        context={
            'site_title':'Plasmid search results',
            'itemlist':object_list,
            'searchcontext':query
            }
    elif type == 'part':
        if query:
            object_list = Magic_pool_part.objects.filter(
                    (Q(name__icontains=query) |
                    Q(description__icontains=query))
                ).distinct().order_by('name')
        else:
            query = ''
            object_list = Magic_pool_part.objects.none()
        template = loader.get_template('magicpool/parts.html')
        context={
            'site_title':'Magic pool part search results',
            'itemlist':object_list,
            'searchcontext':query
            }
    elif type == 'vector':
        if query:
            object_list = Vector.objects.filter(
                    (Q(name__icontains=query) |
                    Q(description__icontains=query))
                ).distinct().order_by('name')
        else:
            query = ''
            object_list = Vector.objects.none()
        template = loader.get_template('magicpool/vectors.html')
        context={
            'site_title':'Vector search results',
            'itemlist':object_list,
            'searchcontext':query
            }
    elif type == 'strain':
        if query:
            object_list = Strain.objects.filter(
                    (Q(name__icontains=query) |
                    Q(amd_number__icontains=query) |
                    Q(species__icontains=query) |
                    Q(description__icontains=query))
                ).distinct().order_by('amd_number')
        else:
            query = ''
            object_list = Strain.objects.none()
        template = loader.get_template('magicpool/strains.html')
        context={
            'site_title':'Strain search results',
            'itemlist':object_list,
            'searchcontext':query
            }
    else:
        return HttpResponse("Not implemented.")
    return HttpResponse(template.render(context, request))

    
