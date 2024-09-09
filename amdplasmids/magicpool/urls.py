from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('search/', views.textsearchform, name='searchform'),
    path('plasmid/<int:plasmid_id>/', views.plasmid_detail, name='plasmiddetails'),
    path('viewer/<int:plasmid_id>/', views.plasmid_viewer, name='plasmidviewer'),
    path('part/<int:part_id>/', views.part_detail, name='partdetails'),
    path('vector/<int:vector_id>/', views.vector_detail, name='vectordetails'),
    path('magicpool/<int:magicpool_id>/', views.magicpool_detail, name='magicpooldetails'),
    path('oligo/<int:oligo_id>/', views.oligo_detail, name='oligodetails'),
    path('strain/<int:strain_id>/', views.strain_detail, name='straindetails'),
    path('plasmids/', views.plasmids, name='plasmidslist'),
    path('parts/', views.parts, name='partslist'),
    path('vectors/', views.vectors, name='vectorslist'),
    path('magicpools/', views.magicpools, name='magicpoolslist'),
    path('strains/', views.strains, name='strainslist'),
    path('nuclsearchform/',views.nucleotidesearchform,name="nucleotidesearchform"),
    path('nuclsearch/',views.nucleotidesearch,name="nucleotidesearch"),
    path('protsearchform/',views.proteinsearchform,name="proteinsearchform"),
    path('protsearch/',views.proteinsearch,name="proteinsearch"),
    path('about/', views.show_help, name='about'),
    path('textsearch/', views.textsearch, name='textsearch'),


    ]
