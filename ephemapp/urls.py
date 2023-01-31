from django.urls import path, include
from . import views

urlpatterns = [
    # ...

    path('sun/', views.multiply_numbers, name='multiply_numbers'),

    # ...
]
