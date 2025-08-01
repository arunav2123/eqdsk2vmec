import unittest
import os
from eqdsk2vmec_converter.eqdsk2vmec import convert_eqdsk_to_vmec

class TestConverter(unittest.TestCase):

    def setUp(self):
        # Create a dummy geqdsk file for testing
        self.test_file = "test.geq"
        with open(self.test_file, "w") as f:
            # This is a mock file, not a real geqdsk file.
            # The freeqdsk library will likely fail to read this.
            # For a real test, a valid geqdsk file is needed.
            f.write("This is a dummy geqdsk file.\n")

    def tearDown(self):
        # Clean up the created files
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
        output_file = f"input.{os.path.splitext(self.test_file)[0]}"
        if os.path.exists(output_file):
            os.remove(output_file)

    def test_conversion(self):
        # This test will likely fail because the dummy file is not a valid geqdsk file.
        # To make this test pass, replace the dummy file with a real geqdsk file.
        try:
            convert_eqdsk_to_vmec(self.test_file)
            # Check if the output file is created
            output_file = f"input.{os.path.splitext(self.test_file)[0]}"
            self.assertTrue(os.path.exists(output_file))
        except Exception as e:
            self.fail(f"Conversion failed with error: {e}")

if __name__ == '__main__':
    unittest.main()


