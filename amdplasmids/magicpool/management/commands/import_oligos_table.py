from django.core.management.base import BaseCommand
from magicpool.util import import_oligos_table

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
            xlsx_path = '/mnt/data/work/Plasmids/datafiles/oligos.xlsx'
        else:
            xlsx_path = options['i']
        import_oligos_table(xlsx_path, options['overwrite'])
