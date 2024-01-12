import importthis as it

# add a record to the database
#it.add_one("Laura","Smith","laura@smith.com")

# will not work if you pass an integer
# instead need to pass as a string
#it.delete_one('6')

# lookup email address record
it.email_lookup("john@codemy.com")


# add many records
'''
stuff = [
	('Brenda', 'Smitherton', 'brenda@smitherton.com'),
	('Joshua', 'Raintree', 'j@rt.com')
	]

it.add_many(stuff)
'''
# show all the records
#it.show_all()




