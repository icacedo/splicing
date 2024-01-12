import sqlite3

conn = sqlite3.connect('customer.db')

c = conn.cursor()

#c.execute("INSERT INTO customers Values ('John', 'Elder', 'john@codemy.com')")

# in order to add another entry, need to re run the script with different insert command

#c.execute("INSERT INTO customers Values ('Tim', 'Smith', 'tim@codemy.com')")

c.execute("INSERT INTO customers Values ('Mary', 'Brown', 'mary@codemy.com')")

# insert many records into the table at the same time

many_customers = [
					('Wes', 'Brown', 'wes@brown.com'), 
					('Steph', 'Kuewa', 'steph@kuewa.com'), 
					('Dan', 'Pas', 'dan@pas.com')
				]

# ? is a place holder
# pass in name of list at the end
c.executemany("INSERT INTO customers VALUES (?,?,?)", many_customers)

conn.commit()

conn.close()
