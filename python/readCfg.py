import re

e = re.compile("(\d+)\s*\n([A-Za-z]{1,2})\s*\n(.*?)(?=\d+\s*\n\s*[A-Za-z]{1,2}\s*\n|\Z)", 
               flags=re.M|re.S)