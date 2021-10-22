import * as React from 'react';
import { ICollection } from './interfaces';

export const CollectionContext = React.createContext<ICollection[]>([]);
