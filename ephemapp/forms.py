from django import forms
import datetime
import pytz

class SextantReadingForm(forms.Form):
    longitude = forms.FloatField(label='Longitude', initial=44.651070)
    latitude = forms.FloatField(label='Latitude', initial=-63.582687)
    date_time = forms.CharField(label='latitude' )
