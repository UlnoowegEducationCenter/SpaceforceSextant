from django.shortcuts import render
from .forms import SextantReadingForm
import ephem

def multiply_numbers(request):
    
    result = "  "

    sun = ephem.Sun()
    halifax = ephem.Observer()
    halifax.lat = "44.651070"
    halifax.lon = "-63.582687"

    halifax.date = "1981/11/10 18:00:00"
    sun.compute(halifax)
    print(sun.dec)

    result = sun.dec

    print("result:", result)
      
              
    return render(request, 'ephemapp/sun.html', { 'result': result})



