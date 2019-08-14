ret_data = {}
# First, we make a list of observation dates from all the bands.
times = [MJD[band] for band in bands]
times = sort(concatenate(times))
# Eliminate repeating days:
gids = concatenate([[1], greater(absolute(times[0:-1] - times[1:]), 0.5)])
times = compress(gids, times)
ret_data['MJD'] = times

# Now loop through the bands and see where we need to fill in data
for band in bands:
	gids = less(absolute(times[:,newaxis] - MJD[band][newaxis,:]), 0.5)
	temp1 = 0.0*times + 99.9
	temp2 = 0.0*times + 99.9
	for i in range(len(gids)):
		if sum(gids[i]) > 0:
			temp1[i] = sum(mags[band]*gids[i])/sum(gids[i])
			temp2[i] =max(sqrt(sum((emags[band]**2)*gids[i])/sum(gids[i])),sqrt(average((temp1[i]-compress(gids[i], mags[band]))**2)))

	ret_data[band] = temp1
	ret_data["e_"+band] = temp2
