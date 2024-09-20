from django.core.management.base import BaseCommand
from magicpool.util import import_strains_table
from amdplasmids.settings import DATA_DIR

class Command(BaseCommand):
    help = '''Imports or updates strain records from xlsx file
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='',
            help='Strain spreadsheet path'
        )
        parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Overwrite existing entries'
        )
    def handle(self, *args, **options):
        if options['i'] == '':
            xlsx_path = os.path.join(DATA_DIR, 'strains.xlsx')
        else:
            xlsx_path = options['i']
        import_strains_table(xlsx_path, options['overwrite'])
