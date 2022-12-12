from django.urls import path
from . import views

app_name = 'vigas'

urlpatterns = [
    path('', views.continua, name='continua'),
    path('calculo-matricial/', views.calcular_matricial, name='calcular_matricial'),
]
