import sys

def tversion(*argv):
    """ Returns load flow version.

    @author: 12345 Canada Inc
    """

    print('arg=',*argv)
    print(len(sys.argv))

    print(sys.argv[0])
    print(len(sys.argv))
    print(str(sys.argv))

    version = {'Name': 'load flow',
           'Version': '0.1.1',
           'Release':  '',
           'Date': '18-Jan-2021'}

    return version


v = tversion('1','2',3)
print (v)
