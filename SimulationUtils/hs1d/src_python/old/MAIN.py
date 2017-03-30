
# coding: utf-8

# In[7]:

import source.ipynb
import TimeUnit.ipynb
import sys

sys.path.append('http://localhost:8888/notebooks')


# In[8]:

class main(object):

    o_Source = Source
    print(o_Source.source(['periodical', None]))

