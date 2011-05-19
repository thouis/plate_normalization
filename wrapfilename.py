import os.path
import re

def wrap_filename(name, postfix):
    base, ext = os.path.splitext(name)
    # if the base ends in _postfix[digits], increment the digits
    m = re.match('(.*)(_%s)_?([0-9]*)'%(re.escape(postfix)), base)
    if m:
        # handle no digits
        if m.group(3):
            count = int(m.group(3)) + 1
        else:
            count = 2
        return m.group(1) + m.group(2) + '_' + str(count) + ext
    return base + '_' + postfix + '_1' + ext



if __name__ == '__main__':
    for name in ['test.xls', 'test_normalized.xls', 'test_normalized_1.xls', 'test_normalized1.xls', 'test', 'test_normalized', 'test_normalized1.xls']:
        print name, wrap_filename(name, 'normalized')
