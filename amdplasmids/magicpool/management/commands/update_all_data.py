from django.core.management.base import BaseCommand
from magicpool.util import update_alldata

class Command(BaseCommand):
    help = '''Update all data from the Google drive
    '''

    def handle(self, *args, **options):
        update_alldata()
