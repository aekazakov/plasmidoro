# Generated by Django 5.0.6 on 2024-07-17 20:55

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('magicpool', '0009_alter_protein_name'),
    ]

    operations = [
        migrations.AddField(
            model_name='plasmid',
            name='magic_pool_part',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='magicpool.magic_pool_part'),
        ),
    ]
