from django.db import models

# Create your models here.
class Contact(models.Model):
    '''
        personal data for contacts
    '''
    name = models.CharField(max_length=255, unique=True)
    email = models.CharField(max_length=255)
    affiliation = models.CharField(max_length=255, blank=True)

    def __str__(self):
        return self.name


class Overhang(models.Model):
    '''
        upstream and downstream overhangs from the magic_pool_design spreadsheet
    '''
    name = models.CharField(max_length=255, unique=True)
    sequence = models.TextField()
    color = models.CharField(max_length=7, default="#FFFFFF")
    
    def __str__(self):
        return self.name
    
class Oligo(models.Model):
    '''
        from the Oligos_and_gBlocks spreadsheet 
    '''
    name = models.CharField(max_length=255, unique=True)
    sequence = models.TextField(blank=True)
    contact = models.ForeignKey(Contact, on_delete=models.SET_NULL, blank=True, null=True)

    def __str__(self):
        return self.name


class Oligo_info(models.Model):
    '''
        unstructured information about oligos
    '''
    oligo = models.ForeignKey(Oligo, on_delete=models.CASCADE)
    param = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    def __str__(self):
        return self.oligo.name + ':' + self.param


class Magic_pool_part_type(models.Model):
    '''
        magic pool part from the magic_pool_design spreadsheet
    '''
    name = models.CharField(max_length=255, unique=True)
    description = models.TextField(blank=True)
    upstream_overhang = models.ForeignKey(Overhang, related_name='upstream_overhang', on_delete=models.SET_NULL, blank=True, null=True)
    downstream_overhang = models.ForeignKey(Overhang, related_name='downstream_overhang', on_delete=models.SET_NULL, blank=True, null=True)

    def __str__(self):
        return self.name


class Magic_pool_part_type_info(models.Model):
    '''
        unstructured information about magic pool parts from the magic_pool_design spreadsheet
    '''
    magic_pool_part_type = models.ForeignKey(Magic_pool_part_type, on_delete=models.CASCADE)
    name = models.CharField(max_length=255)
    description = models.TextField()

    def __str__(self):
        return self.magic_pool_part.name + ': ' +  self.name


class Vector_type(models.Model):
    '''
        Describes design of a magic pool vector
        
    '''
    name = models.CharField(max_length=255, unique=True)
    description = models.TextField(blank=True)

    def __str__(self):
        return self.name


class Vector_type_part(models.Model):
    '''
        ordered list of vector parts
    '''
    vector_type = models.ForeignKey(Vector_type, on_delete=models.CASCADE)
    part_type = models.ForeignKey(Magic_pool_part_type, on_delete=models.CASCADE)
    order = models.PositiveIntegerField()

    class Meta:
        ordering = ['order']
        unique_together = ('vector_type', 'order',)

    def __str__(self):
        return str(self.order) + ': ' + self.part.name


class Vector_type_info(models.Model):
    '''
        unstructured information about a vector
        
    '''
    vector_type = models.ForeignKey(Vector_type, on_delete=models.CASCADE)
    param = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    def __str__(self):
        return self.vector.name + ': ' + self.param


class Feature_type(models.Model):
    '''
        dictionary of feature types
        
    '''
    name = models.CharField(max_length=255, unique=True)

    def __str__(self):
        return self.name


class Drug_marker(models.Model):
    '''
        stores information about drug markers
        
    '''
    name = models.CharField(max_length=255, unique=True)
    drug = models.CharField(max_length=255)
    note = models.CharField(max_length=255, blank=True)

    def __str__(self):
        return self.name


class Strain(models.Model):
    '''
        stores information about strains from the "AMD strain collection" spreadsheet
        
    '''
    amd_number = models.CharField(max_length=255, unique=True)
    description = models.TextField(blank=True)
    plasmid = models.CharField(max_length=255, blank=True)
    species = models.CharField(max_length=255, blank=True)
    name = models.CharField(max_length=255, blank=True)

    def __str__(self):
        return self.amd_number + ' [' + self.species + ']'

        
class Strain_info(models.Model):
    '''
        unstructured information about strains from the "AMD strain collection" spreadsheet
        
    '''
    strain = models.ForeignKey(Strain, on_delete=models.CASCADE)
    param = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    def __str__(self):
        return self.strain.amd_number + ': ' + self.param

        
class Magic_pool(models.Model):
    '''
        magic pool vector from the magic_pool_summary spreadsheet
    '''
    name = models.CharField(max_length=255, unique=True)
    description = models.TextField(blank=True)
    vector_type = models.ForeignKey(Vector_type, on_delete=models.SET_NULL, blank=True, null=True)
    antibiotic_resistance = models.CharField(max_length=255, blank=True)
    strain = models.ForeignKey(Strain, on_delete=models.SET_NULL, blank=True, null=True)

    def __str__(self):
        return self.name


class Plasmid(models.Model):
    '''
        stores information about plasmids from the All_plasmids_Oct19_2023
        
    '''
    name = models.CharField(max_length=255, unique=True)
    amd_number = models.CharField(max_length=255, blank=True)
    description = models.TextField(blank=True)
    sequence = models.TextField(blank=True)
    drug_markers = models.ManyToManyField(Drug_marker)
    magic_pool_designation = models.CharField(max_length=255, blank=True)
    magic_pools = models.ManyToManyField(Magic_pool)
    magic_pool_part = models.ForeignKey(Magic_pool_part_type, on_delete=models.SET_NULL, blank=True, null=True)
    contact = models.ForeignKey(Contact, on_delete=models.SET_NULL, blank=True, null=True)
    sequence_file = models.CharField(max_length=255, blank=True)
    footprint = models.CharField(max_length=32)

    def __str__(self):
        return self.name


class Feature(models.Model):
    '''
        sequence features like genes, IRs, promoters etc.
        
    '''
    name = models.CharField(max_length=255)
    plasmid = models.ForeignKey(Plasmid, on_delete=models.CASCADE)
    feature_type = models.ForeignKey(Feature_type, on_delete=models.SET_NULL, blank=True, null=True)
    sequence = models.TextField(blank=True)
    drug_marker = models.ForeignKey(Drug_marker, on_delete=models.SET_NULL, blank=True, null=True)
    sequence_id = models.CharField(max_length=255)
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    strand = models.IntegerField()
    location_str = models.CharField(max_length=255)
    description = models.TextField(blank=True)

    def __str__(self):
        return self.plasmid.name + ': ' + self.name + ' (' + self.feature_type.name + ')'


class Protein(models.Model):
    '''
        protein sequence
        
    '''
    name = models.CharField(max_length=255)
    sequence = models.TextField()
    function = models.CharField(max_length=255, blank=True)
    feature = models.ForeignKey(Feature, on_delete=models.CASCADE)

    def __str__(self):
        return self.name


class Protein_info(models.Model):
    '''
        protein sequence
        
    '''
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    param = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    def __str__(self):
        return self.protein.name + ': ' + self.param

        
class Plasmid_info(models.Model):
    '''
        unstructured information about plasmids from the All_plasmids_Oct19_2023
        
    '''
    plasmid = models.ForeignKey(Plasmid, on_delete=models.CASCADE)
    param = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    def __str__(self):
        return self.plasmid.name + ': ' + self.param


