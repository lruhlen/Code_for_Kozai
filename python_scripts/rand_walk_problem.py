
def rand_walk(start,left,right,max_count):
    current_pos = start
    count = 0
    
    while (count < max_count) and (current_pos != left) and (current_pos != right):

        go = random.randint(0,2)
        if go == 1:
            current_pos+=1
        else:
            current_pos-=1
        count+=1

        #    print "completed in ",count," steps"
    if current_pos == right:
        return [0,count]
    elif current_pos == left:
        return [1,count]
    else:
        return [-1,count]


def monte_carlo(num_runs):
    start = -12
    left = -50
    right = 10
    max_count = 1e6

    end_left_count = 0
    end_right_count = 0
    successful_run_count = 0
    num_steps_left = []
    num_steps_right = []

    for i in range(num_runs):
        foo = rand_walk(start,left,right,max_count)
        if foo[0] == 0:
            end_right_count+=1
            num_steps_right.append(foo[1])
            successful_run_count+=1
        if foo[0] == 1:
            end_left_count+=1
            num_steps_left.append(foo[1])
            successful_run_count+=1

    print "total of ",successful_run_count," successful runs"
    result = end_left_count * 100.0 / successful_run_count
    print "Ends left first",result,"% of the time"
    clf()
    hist(num_steps_left,alpha=0.4)
    hist(num_steps_right,alpha=0.4)
    return result
