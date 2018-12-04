import unittest

class CommandTest(unittest.TestCase):
    """
    Two ways of running single tests:
    python -m unittest discover -s test -p 'command_test.py'
    python test/command_test.py
    """

    def test_command1(self):
        print("Function 1 is being running")

    def test_command2(self):
        print("Function 2 is being running")

if __name__ == '__main__':
    unittest.main()
