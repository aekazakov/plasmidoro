# Generated by Django 5.0.6 on 2024-06-25 18:49

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('magicpool', '0008_alter_feature_unique_together'),
    ]

    operations = [
        migrations.AlterField(
            model_name='protein',
            name='name',
            field=models.CharField(max_length=255),
        ),
    ]