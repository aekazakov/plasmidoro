from django.core.management.base import BaseCommand
from magicpool.util import import_oligos_table
from amdplasmids.settings import DATA_DIR

class Command(BaseCommand):
    help = '''Imports or updates oligo records from xlsx file
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='',
            help='Oligos spreadsheet path'
        )
        parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Overwrite existing entries'
        )
    def handle(self, *args, **options):
        if options['i'] == '':
            xlsx_path = os.path.join(DATA_DIR, 'oligos.xlsx')
        else:
            xlsx_path = options['i']
        import_oligos_table(xlsx_path, options['overwrite'])
