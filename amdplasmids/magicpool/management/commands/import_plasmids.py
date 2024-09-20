from django.core.management.base import BaseCommand
from magicpool.util import import_plasmids
from amdplasmids.settings import DATA_DIR

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
            work_dir = os.path.join(DATA_DIR, 'plasmid_maps')
        else:
            work_dir = options['i']
        print(work_dir)
        import_plasmids(work_dir, options['overwrite'])
