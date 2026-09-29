class Dragger:
	def __init__(self, fig, release_callback=None, move_callback=None):
		self.fig = fig
		self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
		self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
		self.fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
		self.release_callback = release_callback
		self.move_callback = move_callback

		self.mouse_held = False

	def button_press_callback(self, event):
		for artist in self.fig.axes[0].get_xticklabels():
			if not artist.contains(event)[0]:
				continue
			self.mouse_held = True
			self.selected_label = artist
			self.selected_label_position_old = artist.get_position()[0]
			self.label_position_limits = self.fig.axes[0].get_xlim()
			self.xtick_positions = {}			
			self.xtick_labels = {}

			xticks = self.fig.axes[0].get_xticklabels()
			for xtick in xticks:
				if xtick == self.selected_label:
					self.xtick_positions[xtick] = self.selected_label_position_old
				else:
					self.xtick_positions[xtick] = xtick.get_position()[0]
				self.xtick_labels[xtick] = xtick.get_text()

			self.xtick_position_diff = list(sorted(self.xtick_positions.values()))[1] - list(sorted(self.xtick_positions.values()))[0]
			return

	def button_release_callback(self, event):
		if self.mouse_held:
			self.xtick_positions[self.selected_label] = self.selected_label_position_old
			self.update_xtick_positions()
			self.update_fig()

			self.mouse_held = False
			self.selected_label = None
			self.selected_label_position = None
			self.label_position_limits = None

			if self.release_callback:
				self.release_callback(self.xtick_positions)

	def motion_notify_callback(self, event):
		if not self.mouse_held:
			return

		# we calculate the visible position of the xtick label
		new_x = self.fig.axes[0].transData.inverted().transform((event.x, 0))[0]
		if new_x > self.label_position_limits[1]:
			new_x = self.label_position_limits[1]
		if new_x < self.label_position_limits[0]:
			new_x = self.label_position_limits[0]

		# set it here temporarily. We will reinstate the old position when the mouse is released
		self.xtick_positions[self.selected_label] = new_x

		# now we also calculate the distances to all other ones
		# and the distance to the original position
		sdistances = {'old': new_x - self.selected_label_position_old}
		for xtick, x in self.xtick_positions.items():
			if xtick == self.selected_label:
				continue
			sdistances[xtick] = new_x - x

		# find which item is closest to the currently moving item
		closest_item, lowest_distance = min(sdistances.items(), key=lambda k: abs(k[1]))
		# and obtain which items are inbetween the old position and the moving item
		inbetween_items = []
		if sdistances['old'] > 0:
			# define a range within which the item distances should be
			# (0, pos. number)
			rang = (0, sdistances['old'])
			# the shift tells you which way to move the items if they are between the
			# old and current positions
			shift = -self.xtick_position_diff
		else:
			# (neg. number, 0)
			rang = (sdistances['old'], 0)
			shift = self.xtick_position_diff

		# check now which items are between the old an current position
		for xtick, sdist in sdistances.items():
			if xtick == 'old':
				continue

			# also add the closest item, which may or may not be inbetween
			if xtick == closest_item:
				inbetween_items.append(xtick)
				continue

			if rang[0] < sdist < rang[1]:
				inbetween_items.append(xtick)

		# if we do not have any inbetween items we skip this part
		if len(inbetween_items) > 0:

			# if we are close to another item we set the current item's position to its
			self.selected_label_position_old = self.xtick_positions[closest_item]

			# and then shift all inbetween items
			for inbetween_item in inbetween_items:
				self.xtick_positions[inbetween_item] += shift

		self.update_xtick_positions()
		self.update_fig()

		if self.move_callback:
			self.move_callback(self.xtick_positions)

	def update_xtick_positions(self):
		pos = list(self.xtick_positions.values())
		labs = list(self.xtick_labels.values())
		self.fig.axes[0].set_xticks(pos, labs)

	def update_fig(self):
		self.fig.canvas.draw_idle()

