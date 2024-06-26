from django.core.management.base import BaseCommand
from magicpool.util import import_plasmids_table

class Command(BaseCommand):
    help = '''Imports or updates plasmid records from xlsx file
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='',
            help='Top directory with plasmid files'
        )
        parser.add_argument(
            '--overwrite',
            action='store_true',
            help='Overwrite existing entries'
        )
    def handle(self, *args, **options):
        if options['i'] == '':
            xlsx_path = '/mnt/data/work/Plasmids/datafiles/plasmid_maps/All_plasmids_Oct19_2023.xlsx'
        else:
            xlsx_path = options['i']
        import_plasmids_table(xlsx_path)
