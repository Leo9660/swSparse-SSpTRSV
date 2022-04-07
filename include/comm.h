#ifndef STRUCT_COMM_
#define STRUCT_COMM_

//Data type
#define TYPE double
#define TYPE_V doublev8
#define MEM_FENCE asm volatile("memb":::"memory")

//communication
static inline void* remote_ldm_addr(volatile void* p, uintptr_t id) {
  uintptr_t x = (uintptr_t)p, y = (uintptr_t)id << 20, z = 1ull << 45;
  return (void*)(x | y | z);
}

static inline void rstore8(TYPE* buf, volatile void* p, uintptr_t id)
{
    TYPE_V vbuf;
    TYPE *addr = remote_ldm_addr(p, id);
    simd_loadu(vbuf, buf);
    simd_storeu(vbuf, addr);
}

#endif