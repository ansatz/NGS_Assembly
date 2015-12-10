import functools
import inspect

#http://stackoverflow.com/questions/9906144/python-decorate-a-class-by-defining-the-decorator-as-a-class
#http://stackoverflow.com/questions/739654/how-can-i-make-a-chain-of-function-decorators-in-python
class object_proxy(object): 

    def __init__(self, wrapped):
        self.wrapped = wrapped
        try:
            self.__name__= wrapped.__name__
        except AttributeError:
            pass 

    @property
    def __class__(self):
        return self.wrapped.__class__  

    def __getattr__(self, name):
        return getattr(self.wrapped, name)


class bound_function_wrapper(object_proxy): 

    def __init__(self, wrapped, instance, wrapper, binding, parent):
        super(bound_function_wrapper, self).__init__(wrapped)
        self.instance = instance
        self.wrapper = wrapper
        self.binding = binding
        self.parent = parent 

    def __call__(self, *args, **kwargs):
        if self.binding == 'function':
            if self.instance is None:
                instance, args = args[0], args[1:]
                wrapped = functools.partial(self.wrapped, instance)
                return self.wrapper(wrapped, instance, args, kwargs)
            else:
                return self.wrapper(self.wrapped, self.instance, args, kwargs)
        else:
            instance = getattr(self.wrapped, '__self__', None)
            return self.wrapper(self.wrapped, instance, args, kwargs) 

    def __get__(self, instance, owner):
        if self.instance is None and self.binding == 'function':
            descriptor = self.parent.wrapped.__get__(instance, owner)
            return bound_function_wrapper(descriptor, instance, self.wrapper,
                    self.binding, self.parent)
        return self 

class function_wrapper(object_proxy): 

    def __init__(self, wrapped, wrapper):
        super(function_wrapper, self).__init__(wrapped)
        self.wrapper = wrapper
        if isinstance(wrapped, classmethod):
            self.binding = 'classmethod'
        elif isinstance(wrapped, staticmethod):
            self.binding = 'staticmethod'
        else:
            self.binding = 'function' 

    def __get__(self, instance, owner):
        wrapped = self.wrapped.__get__(instance, owner)
        return bound_function_wrapper(wrapped, instance, self.wrapper,
                self.binding, self) 

    def __call__(self, *args, **kwargs):
        return self.wrapper(self.wrapped, None, args, kwargs)

#decorator factory
def decorator(wrapper):
    @functools.wraps(wrapper)
    def _decorator(wrapped):
        return function_wrapper(wrapped, wrapper)
    return _decorator

@decorator
def universal(wrapped, instance, args, kwargs):
    if instance is None:
        if inspect.isclass(wrapped):
            # Decorator was applied to a class.
            return wrapped(*args, **kwargs)
        else:
            # Decorator was applied to a function or staticmethod.
            return wrapped(*args, **kwargs)
    else:
        if inspect.isclass(instance):
            # Decorator was applied to a classmethod.
            return wrapped(*args, **kwargs)
        else:
            # Decorator was applied to an instancemethod.
            return wrapped(*args, **kwargs)

@decorator
def my_function_wrapper(wrapped, instance, args, kwargs):
    print('WRAPPED', wrapped)
    print('INSTANCE', instance)
    print('ARGS', args)
    return wrapped(*args, **kwargs) 


@universal
@my_function_wrapper
def function(a,b):
    pass

