from django import forms
class ObservationForm(forms.Form):
    star_number = forms.IntegerField(min_value=0)
    elevation = forms.FloatField()
