from django.core.management.base import BaseCommand
from magicpool.util import import_magic_pool_types
from amdplasmids.settings import DATA_DIR


class Command(BaseCommand):
    help = '''Imports or updates magic pool parts and vector designs from an xlsx file
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='',
            help='xlsx file path'
        )
    def handle(self, *args, **options):
        if options['i'] == '':
            xlsx_path = os.path.join(DATA_DIR, 'plasmid_maps', 'magicpool_vector_designs', 'magic_pool_design.xlsx')
        else:
            xlsx_path = options['i']
        import_magic_pool_types(xlsx_path)
