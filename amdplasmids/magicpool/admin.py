from django.contrib import admin

# Register your models here.
from magicpool.models import *

admin.site.register(Contact)
admin.site.register(Overhang)
admin.site.register(Oligo)
admin.site.register(Oligo_info)
admin.site.register(Magic_pool)
admin.site.register(Magic_pool_part_type)
admin.site.register(Magic_pool_part_type_info)
admin.site.register(Vector_type_part)
admin.site.register(Vector_type)
admin.site.register(Feature_type)
admin.site.register(Protein)
admin.site.register(Protein_info)
admin.site.register(Drug_marker)
admin.site.register(Plasmid)
admin.site.register(Feature)
admin.site.register(Plasmid_info)
admin.site.register(Strain)
admin.site.register(Strain_info)
