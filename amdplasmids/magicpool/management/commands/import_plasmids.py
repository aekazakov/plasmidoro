from django.core.management.base import BaseCommand
from magicpool.util import update_plasmids

class Command(BaseCommand):
    help = '''Imports or updates plasmid records from directory tree
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
            work_dir = '/mnt/data/work/Plasmids/datafiles/plasmid_maps'
        else:
            work_dir = options['i']
        print(work_dir)
        update_plasmids(work_dir, options['overwrite'])
