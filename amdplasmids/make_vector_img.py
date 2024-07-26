"""
    Generate vector design image with plasmidcanvas
    
    Invoke: exec(open('make_vector_img.py').read())
"""
from plasmidcanvas.plasmid import *
from plasmidcanvas.feature import *
from magicpool.models import Vector
from magicpool.models import Vector_part
from magicpool.models import Plasmid as Plasmidoro_plasmid
from magicpool.models import Plasmid_info
from magicpool.models import Feature
from magicpool.models import Vector

def process_plasmid_no_overhangs(plasmid):
	plasmid_info = Plasmid_info.objects.filter(plasmid=plasmid, param='Sequence from BbsI sites')
	if plasmid_info.exists():
		sequence = str(plasmid_info[0].value)
		sequence = sequence.lower()
		upstream_overhang = plasmid.magic_pool_part.upstream_overhang.sequence.lower()
		downstream_overhang = plasmid.magic_pool_part.downstream_overhang.sequence.lower()
		upstream_overhang_pos = sequence.find(upstream_overhang)
		if upstream_overhang_pos == -1:
			raise Exception('ERROR: upstream overhang not found in  ' + plasmid.name)
		downstream_overhang_pos = sequence.rfind(downstream_overhang)
		if downstream_overhang_pos == -1:
			raise Exception('ERROR: downstream overhang not found in  ' + plasmid.name)
		return [], downstream_overhang_pos - upstream_overhang_pos - len(upstream_overhang) + 1
	else:
		raise Exception('ERROR: overhang not found in  ' + plasmid.name)
		

def get_vector_part_features(plasmid, overhang_size):
	min_plasmid_overhang_pos = None
	max_plasmid_overhang_pos = 0
	vector_part_size = 0
	for feature in Feature.objects.filter(plasmid=plasmid).all():
		if feature.name.startswith('4b'):
			if min_plasmid_overhang_pos is None:
				# initialize vars
				min_plasmid_overhang_pos = feature.start
				max_plasmid_overhang_pos = feature.start
			elif min_plasmid_overhang_pos > feature.start:
				min_plasmid_overhang_pos = feature.start
			elif max_plasmid_overhang_pos < feature.start:
				max_plasmid_overhang_pos = feature.start
	if min_plasmid_overhang_pos is None:
		return process_plasmid_no_overhangs(plasmid)
	if min_plasmid_overhang_pos == max_plasmid_overhang_pos:
		return process_plasmid_no_overhangs(plasmid)
	# Choose inward or outward direction
	vector_part_size = max_plasmid_overhang_pos - min_plasmid_overhang_pos + 1 - overhang_size
	vector_features = []
	shift = min_plasmid_overhang_pos + overhang_size - 1
	for feature in Feature.objects.filter(plasmid=plasmid).all():
		if feature.name.startswith('4b'):
			continue
		if feature.name.startswith('BsmBI'):
			continue
		if feature.name.startswith('changed base'):
			continue
		if feature.name.startswith('base change'):
			continue
		if feature.start > shift and feature.end < max_plasmid_overhang_pos:
			if feature.strand == 1:
				vector_features.append([feature.name, feature.start - shift, feature.end - shift])
			else:
				vector_features.append([feature.name, feature.end - shift, feature.start - shift])
	return vector_features, vector_part_size


#vector_name = 'BarTn7_onepot'

for vector_name in Vector.objects.values_list('name', flat=True):
	print('    PROCESSING ' + vector_name)
	vector_parts = []
	vector_features = []
	overhang_size = 4
	vector_size = overhang_size

	# Collect the data
	vector = Vector.objects.get(name=vector_name)
	skip_vector = False
	for vector_part in Vector_part.objects.filter(vector=vector):
		magic_pool_part = vector_part.part
		if Plasmidoro_plasmid.objects.filter(magic_pool_part=magic_pool_part).exists():
			plasmid = Plasmidoro_plasmid.objects.filter(magic_pool_part=magic_pool_part)[0]
			vector_parts.append([magic_pool_part.upstream_overhang.name + ' overhang', vector_size - overhang_size + 1, vector_size])
			features, part_size = get_vector_part_features(plasmid, overhang_size)
			vector_parts.append([magic_pool_part.name, vector_size + 1, vector_size + part_size])
			for feature in features:
				vector_features.append([feature[0], feature[1] + vector_size, feature[2] + vector_size])
			vector_size += part_size + overhang_size
		else:
			skip_vector = True
			print('ERROR: vector part for ' + magic_pool_part.name + ' not found')
	if skip_vector:
		continue
	print(vector_parts)
	print(vector_features)

	COLORS = ['red', 'green', 'blue', 'orange', 'darkred', 'purple']
	# Now let's draw!
	plasmid_img = Plasmid(vector_name.replace('_', ''), vector_size)
	plasmid_img.set_feature_label_font_size(6)
	color_ind = 0
	for vector_part in vector_parts:
		if vector_part[0].endswith('overhang'):
			label = SinglePairLabel(vector_part[0], vector_part[1])
			label.set_font_color("black")
			label.set_font_size(4)
			plasmid_img.add_feature(label)
		else:
			some_part = RectangleFeature(vector_part[0], vector_part[1], vector_part[2])
			if vector_part[2] - vector_part[1] > 1000:
				some_part.set_label_styles(["on-circle"])
			some_part.set_color(COLORS[color_ind])
			plasmid_img.add_feature(some_part)
			color_ind += 1
			if color_ind > len(COLORS) - 1:
				color_ind = 0
		
	for vector_feature in vector_features:
		some_feature = ArrowFeature(vector_feature[0], vector_feature[1], vector_feature[2], direction=1)
		plasmid_img.add_feature(some_feature)
	plasmid_img.save_to_file(vector_name + ".svg")
