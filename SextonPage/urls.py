from django.urls import path
from ephemapp import views
from django.contrib import admin

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.multiply_numbers, name='multiply_numbers'),
    path('sun/', views.multiply_numbers, name='multiply_numbers'),
]
