from django import forms
import datetime
import pytz

class SextantReadingForm(forms.Form):
    longitude = forms.CharField(required= 'False',label='longitude')
    latitude = forms.CharField(label='latitude')
    date_time = forms.CharField(label='date_time')
