from functools import wraps
import time


def timer(func):
    '''
    Using decorators to do the db cleaning ...
     - This deletes any extant db before the test
     - And deletes any extant db after the test
    '''
    @wraps(func)
    def inner_function(*args, **kwargs):
    
        # Time before
        start_time = time.clock()
            
        # Run target function
        result = func(*args, **kwargs)

        # Time after (and hence delta-t)
        end_time = time.clock()
        delta_t = end_time - start_time
        
        # print / log
        print(  f'\n',
                f'[func=      ]{func.__name__}\n',
                f'[start_time=]{start_time}\n',
                f'[end_time=  ]{end_time}\n',
                f'[delta_t=   ]{delta_t}\n',
            )
        
        # Return any result
        return result

    return inner_function

