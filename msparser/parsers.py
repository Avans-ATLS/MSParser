import re

class Formula:
    def __init__(self, formula: str) -> None:
        """Read a chemical formula and parse it into elements and counts

        Args:
            formula (str): chemical formula
        """
        # Dictionary to store elements and their counts, isotopes handled separately
        self.elements = {}
        # Regular expression pattern for elements and isotopic notation (optional)
        regexpattern = r"(\[?\d*[A-Z][a-z]*\]?)(\d*)"
        matches = re.findall(regexpattern, formula)
        
        for element, count in matches:
            count = int(count) if count else 1
            
            # Check if it’s an isotope (if it contains square brackets, e.g., [13C])
            if element.startswith('['):
                # Separate the isotope notation and base element, e.g., '[13C]' becomes '13' and 'C'
                isotope = re.findall(r"\[(\d+)([A-Z][a-z]*)\]", element)[0]
                base_element = isotope[1]
                isotope_number = int(isotope[0])
                # Add isotope count to the element in a tuple
                if base_element not in self.elements:
                    self.elements[base_element] = {'standard': 0, 'isotopes': {}}
                if isotope_number in self.elements[base_element]['isotopes']:
                    self.elements[base_element]['isotopes'][isotope_number] += count
                else:
                    self.elements[base_element]['isotopes'][isotope_number] = count
            else:
                # Handle standard elements
                if element not in self.elements:
                    self.elements[element] = {'standard': 0, 'isotopes': {}}
                self.elements[element]['standard'] += count
    
    def change_element_count(self, element: str, modifier: str, count: int, isotope: int = None) -> None:
        """Change the count of an element or isotope in the formula

        Args:
            element (str): The element to change
            modifier (str): add ('+') or subtract ('-')
            count (int): count of elements or isotopes to add or subtract
            isotope (int, optional): If modifying an isotope, provide the isotope number. Defaults to None.

        Raises:
            ValueError: If invalid modifier is provided
        """
        if isotope:
            # Modify isotopic count
            if element not in self.elements:
                self.elements[element] = {'standard': 0, 'isotopes': {}}
            if modifier == "+":
                self.elements[element]['isotopes'][isotope] = self.elements[element]['isotopes'].get(isotope, 0) + count
            elif modifier == "-":
                self.elements[element]['isotopes'][isotope] -= count
                if self.elements[element]['isotopes'][isotope] <= 0:
                    del self.elements[element]['isotopes'][isotope]  # Remove isotope if count goes to zero
            else:
                raise ValueError("Modifier must be '+' or '-'")
        else:
            # Modify standard element count
            if modifier == "+":
                self.elements[element]['standard'] += count
            elif modifier == "-":
                self.elements[element]['standard'] -= count
                if self.elements[element]['standard'] <= 0:
                    del self.elements[element]  # Remove element if count goes to zero
            else:
                raise ValueError("Modifier must be '+' or '-'")
    
    def to_string(self):
        """Rebuilds the chemical formula as a string, ensuring CHNO are in order,
        followed by other elements alphabetically, and isotopes are placed after standard elements"""
        
        formula = ""
        
        # Custom element order: CHNO, then others alphabetically
        custom_order = ['C', 'H', 'N', 'O']
        all_elements = list(self.elements.keys())
        
        # Sort elements by custom order (CHNO first), then alphabetically for the rest
        sorted_elements = sorted(all_elements, key=lambda e: (e not in custom_order, e))
        
        for element in sorted_elements:
            data = self.elements[element]
            
            # First, add the standard element if it has a count > 0
            if data['standard'] > 0:
                formula += element + (str(data['standard']) if data['standard'] > 1 else "")
            
            # Add isotopes, ensuring they come after the standard form
            for isotope, count in sorted(data['isotopes'].items()):
                formula += f"[{isotope}{element}]{count if count > 1 else ''}"
        
        return formula
