function progress_bar(progress::Number, width=displaysize(stdout)[2], keep_bar=true)
	progress = float(progress)
	if progress < 0
		progress = 0.0
	elseif progress > 1
		progress = 1.0
	end
	width = displaysize(stdout)[2]
	if progress == 1
		if keep_bar
			println("\r" * "█"^round(Int, progress*width) * "━"^round(Int, (1-progress)*width))
		else
			print("\r" * "█"^round(Int, progress*width) * "━"^round(Int, (1-progress)*width))
		end
	else
		print("\r" * "█"^round(Int, progress*width) * "━"^round(Int, (1-progress)*width))
	end
end

function simulate_bar(secs)
	t = 0
	dt = 0.01
	while t <= secs
		sleep(dt)
		progress_bar(t/secs)
		t += dt
	end
	progress_bar(1)
end
simulate_bar(1.0)
