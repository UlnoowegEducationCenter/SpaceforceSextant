from django.shortcuts import render
from .forms import SextantReadingForm
import ephem

def multiply_numbers(request):
    
    result = " "

    if request.method == 'POST':
        form = SextantReadingForm(request.POST)
        if form.is_valid():
            date_time = form.cleaned_data['date_time']
            longitude = form.cleaned_data['longitude']
            latitude = form.cleaned_data['latitude']
     

            try:
                sun = ephem.Sun()
                halifax = ephem.Observer()
              

                halifax.date = date_time
                sun.compute(halifax)
                print(sun.dec)

                result = sun.dec

                print("result:", result)
                
            except Exception as e:
                print("Error:", e)
    else:
        form = SextantReadingForm()
    return render(request, 'ephemapp/sun.html', {'form': form, 'result': result})



