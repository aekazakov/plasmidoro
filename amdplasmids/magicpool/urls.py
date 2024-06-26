from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('search/', views.textsearchform, name='searchform'),
    path('plasmid/<int:plasmid_id>/', views.plasmid_detail, name='plasmiddetails'),
    path('viewer/<int:plasmid_id>/', views.plasmid_viewer, name='plasmidviewer'),
    path('part/<int:part_id>/', views.part_detail, name='partdetails'),
    path('vector/<int:vector_id>/', views.vector_detail, name='vectordetails'),
    path('plasmids/', views.plasmids, name='plasmidslist'),
    path('parts/', views.parts, name='partslist'),
    path('vectors/', views.vectors, name='vectorslist'),
    path('nuclsearchform/',views.nucleotidesearchform,name="nucleotidesearchform"),
    path('nuclsearch/',views.nucleotidesearch,name="nucleotidesearch"),
    path('protsearchform/',views.proteinsearchform,name="proteinsearchform"),
    path('protsearch/',views.proteinsearch,name="proteinsearch"),
    path('help/', views.show_help, name='help'),
    path('textsearch/', views.textsearch, name='textsearch'),


    ]
