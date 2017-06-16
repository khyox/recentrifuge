"""mypy test: location for undefined references in f-strings"""

from typing import List

my_list: List[str] = ['Anna', 'Jonathan', 'Sandra', 'Joseph']
print(f'My list has {len(my_list)} elements.')

print('The list has %i elements.' % len(my_wrong_list))
print(f'The list has {len(another_wrong_list)} elements.')
print(f'The list has {missing_integer} elements.')
print('The list has {0:d} elements.'.format(len(yet_another_wrong_list)))
