from django.contrib import admin
from django.urls import path
from sextantCalculate import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('calculate_position/', views.calculate_position_view, name='calculate_position'),
]
