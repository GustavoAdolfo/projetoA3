#!/bin/bash
rm -fr .vscode
rm -fr ENV_DIR
rm db.sqlite3 db.sqlite3.old geckodriver.log requirements-dev.txt
# rm calc_vigas/settings.py
# mv calc_vigas/settings.py.prd calc_vigas/settings.py
# rm calc_vigas/.env
# mv calc_vigas/.env.prd calc_vigas/.env
# python ./createdb.py
python manage.py makemigrations
python manage.py migrate
python manage.py compilescss
python manage.py collectstatic --no-input
echo "from django.contrib.auth import get_user_model; User = get_user_model(); User.objects.create_superuser(email='calc_vigas@livrosviajantes.com.br', password='S3nh@root')" | python manage.py shell;

# gunicorn --bind 0.0.0.0:80 --workers 3 calc_vigas.wsgi:application --log-level info

python manage.py runserver 0.0.0.0:80

exec "$@"