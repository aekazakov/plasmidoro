from django.core.management.base import BaseCommand
from django.core.mail import mail_admins
from magicpool.util import update_alldata

class Command(BaseCommand):
    help = '''Update all data from the Google drive
    '''
    def add_arguments(self, parser):
        parser.add_argument('--replace', action='store_true', help='Update existing plasmids if plasmid map files has changed')

    def handle(self, *args, **options):
        report = update_alldata(options['replace'])
        mail_admins('Plasmidoro: update all data finished', '\n'.join(report))