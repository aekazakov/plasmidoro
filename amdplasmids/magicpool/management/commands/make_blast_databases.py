from django.core.management.base import BaseCommand
from magicpool.util import make_blast_databases

class Command(BaseCommand):
    help = '''Generates blast databases
    '''

    def handle(self, *args, **options):
        make_blast_databases()
