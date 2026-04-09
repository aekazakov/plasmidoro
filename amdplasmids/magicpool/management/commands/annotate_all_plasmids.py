from django.core.management.base import BaseCommand
from magicpool.plannotate import application
from magicpool.models import Plasmid

class Command(BaseCommand):
    help = '''Annotates all plasmids
    '''
    def handle(self, *args, **options):
        plasmids = Plasmid.objects.all()
        application(plasmids)
