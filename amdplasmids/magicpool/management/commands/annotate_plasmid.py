from django.core.management.base import BaseCommand
from magicpool.plannotate import application
from magicpool.models import Plasmid

class Command(BaseCommand):
    help = '''Annotates a plasmid
    '''

    def add_arguments(self, parser):
        parser.add_argument(
            '-i',
            default='',
            help='Plasmid name'
        )
    def handle(self, *args, **options):
        if options['i'] != '':
            plasmid_name = options['i']
            plasmids = Plasmid.objects.filter(name = plasmid_name)
            application(plasmids)
        else:
            print('Plasmid name should not be empty')
