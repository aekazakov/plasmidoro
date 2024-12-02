from django.core.management.base import BaseCommand
from magicpool.util import update_alldata

class Command(BaseCommand):
    help = '''Update all data from the Google drive
    '''
    def add_arguments(self, parser):
        parser.add_argument('--replace', action='store_true', help='Update existing plasmids if plasmid map files has changed')

    def handle(self, *args, **options):
        update_alldata(options['replace'])
