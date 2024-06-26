from django.shortcuts import render
from django.template import loader
from django.http import HttpResponse
from django.db.models import Q
from magicpool.models import *

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
        'info':Magic_pool_part_info.objects.filter(magic_pool_part=part_id)
        }
    return HttpResponse(template.render(context, request))


def vector_detail(request, vector_id):
    '''
        Displays vector design page
    '''
    template = loader.get_template('magicpool/vector.html')
    vector = Vector.objects.get(id=vector_id)
    context = {'vector':vector,
        'info':Vector_info.objects.filter(vector=vector_id),
        'vector_parts':Vector_part.objects.filter(vector=vector_id).select_related('part')
        }
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


def nucleotidesearchform(request):
    '''
        Displays nucleotide sequence search form
    '''
    return render(request,'magicpool/nucleotidesearchform.html')


def nucleotidesearch(request):
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
        Displays protein sequence search form
    '''
    #return render(request,'magicpool/proteinsearchform.html')
    return HttpResponse("Not implemented.")


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
    else:
        return HttpResponse("Not implemented.")
    return HttpResponse(template.render(context, request))

    
