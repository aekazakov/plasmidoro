from django.core.management.base import BaseCommand
from magicpool.util import import_magic_pool

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
            xlsx_path = '/mnt/data/work/Plasmids/datafiles/plasmid_maps/magicpool_vector_designs/magic_pool_design.xlsx'
        else:
            xlsx_path = options['i']
        import_magic_pool(xlsx_path)
